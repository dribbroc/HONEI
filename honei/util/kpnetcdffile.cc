/* vim: set sw=4 sts=4 et foldmethod=syntax : */
/*
 * KPNetCDFFile.cpp
 *
 *  Created on: Apr 9, 2010
 *      Author: babrodtk
 */

#include <honei/util/kpnetcdffile.hh>
#include <honei/util/attributes.hh>

#include <iostream>
#include <cstdlib>

using namespace std;
using namespace tr1;

//Open a new file?
KPNetCDFFile::KPNetCDFFile(std::string filename, const KPInitialConditions& init) {
    cout << "Opening new file '" << filename << "' " << endl;
#ifdef USE_NETCDF_COMPRESSION
    file.reset(new NcFile(filename.c_str(), NcFile::New, 0, 0, NcFile::Netcdf4));
#else
    file.reset(new NcFile(filename.c_str(), NcFile::New));
#endif
    if (!file->is_valid()) {
        cout << "Could not create '" << filename << "'." << endl;
        cout << "Check that it does not exist, or that your disk is full." << endl;
        exit(-1);
    }
    writeDims(init);
    writeVars(init);
    writeAtts(init);
    nt = 0;
    time_offset = 0;
    writeInitialConditions(init);
}

//Or try opening existing file
KPNetCDFFile::KPNetCDFFile(std::string filename, KPInitialConditions& init,
        shared_ptr<Field>& init_B,
        shared_ptr<Field>& init_U1,
        shared_ptr<Field>& init_U2,
        shared_ptr<Field>& init_U3) {
    cout << "Opening existing file '" << filename << "' " << endl;
    file.reset(new NcFile(filename.c_str(), NcFile::Write));
    if (!file->is_valid())
        file.reset(new NcFile(filename.c_str(), NcFile::Write, 0, 0, NcFile::Netcdf4));
    if (!file->is_valid()) {
        cout << "Could not open..." << endl;
        exit(-1);
    }
    readDims(init);
    readVars(init);
    readAtts(init);
    nt = layout.dims.t->size() - 1;
    layout.vars.t->set_cur(nt);
    layout.vars.t->get(&time_offset, 1);
    readInitialConditions(init, init_B, init_U1, init_U2, init_U3);
}

KPNetCDFFile::~KPNetCDFFile() {
    file->sync();
    if (!file->close()) {
        cout << "Error: Couldn't close NetCDF file!" << endl;
    }
}

void KPNetCDFFile::writeDims(const KPInitialConditions& init) {
    float* tmp = new float[max(init.getBNx(), init.getBNy())];

    //Create dimensions
    layout.dims.i = file->add_dim("I", init.getBNx());
    layout.dims.j = file->add_dim("J", init.getBNy());
    layout.dims.x = file->add_dim("X", init.getNx());
    layout.dims.y = file->add_dim("Y", init.getNy());
    layout.dims.t = file->add_dim("T");

    //Create indexing variables
    layout.vars.i = file->add_var("I", ncFloat, layout.dims.i);
    layout.vars.j = file->add_var("J", ncFloat, layout.dims.j);
    layout.vars.i->add_att("description", "Longitudal coordinate for values given at grid cell intersections");
    layout.vars.j->add_att("description", "Latitudal coordinate for values given at grid cell intersections");

    layout.vars.x = file->add_var("X", ncFloat, layout.dims.x);
    layout.vars.y = file->add_var("Y", ncFloat, layout.dims.y);
    layout.vars.x->add_att("description", "Longitudal coordinate for values given at grid cell centers");
    layout.vars.y->add_att("description", "Latitudal coordinate for values given at grid cell centers");

    layout.vars.t = file->add_var("T", ncFloat, layout.dims.t);
    layout.vars.t->add_att("description", "Time");

    //Write contents of spatial variables
    for (uint i=0; i<init.getBNx(); ++i)
        tmp[i] = i * init.getDx();
    layout.vars.i->put(tmp, init.getBNx());

    for (uint i=0; i<init.getBNy(); ++i)
        tmp[i] = i * init.getDy();
    layout.vars.j->put(tmp, init.getBNy());

    for (uint i=0; i<init.getNx(); ++i)
        tmp[i] = (i+0.5) * init.getDx();
    layout.vars.x->put(tmp, init.getNx());

    for (uint i=0; i<init.getNy(); ++i)
        tmp[i] = (i+0.5) * init.getDy();
    layout.vars.y->put(tmp, init.getNy());

    delete [] tmp;
    file->sync();
}

void KPNetCDFFile::readDims(HONEI_UNUSED KPInitialConditions& init) {
    //Get dimensions
    layout.dims.i = file->get_dim("I");
    layout.dims.j = file->get_dim("J");
    layout.dims.x = file->get_dim("X");
    layout.dims.y = file->get_dim("Y");
    layout.dims.t = file->get_dim("T");

    //Get indexing variables
    layout.vars.i = file->get_var("I");
    layout.vars.j = file->get_var("J");

    layout.vars.x = file->get_var("X");
    layout.vars.y = file->get_var("Y");

    layout.vars.t = file->get_var("T");
}

void KPNetCDFFile::writeVars(const KPInitialConditions& init) {
    //Create initial condidion variables
    if (init.getBLocation() == KPInitialConditions::GRID_CELL_INTERSECTION) {
        layout.vars.init_B = file->add_var("init_B", ncFloat, layout.dims.j, layout.dims.i);
        layout.vars.init_B->add_att("description", "Initial conditions for bathymetry");
    }
    else {
        cout << "Error: Go fix it you fool!" << endl;
        exit(-1);
    }

    //Create the timestep variables
    layout.vars.U1 = file->add_var("U1", ncFloat, layout.dims.t, layout.dims.y, layout.dims.x);
    layout.vars.U2 = file->add_var("U2", ncFloat, layout.dims.t, layout.dims.y, layout.dims.x);
    layout.vars.U3 = file->add_var("U3", ncFloat, layout.dims.t, layout.dims.y, layout.dims.x);

    layout.vars.U1->add_att("description", "Water elevation");
    layout.vars.U2->add_att("description", "Longitudal water discharge");
    layout.vars.U3->add_att("description", "Latitudal water discharge");


#ifdef USE_NETCDF_COMPRESSION
    nc_def_var_deflate(file->id(), layout.B_init->id(), 1, 1, 2);
    nc_def_var_deflate(file->id(), layout.U1_init->id(), 1, 1, 2);
    nc_def_var_deflate(file->id(), layout.U2_init->id(), 1, 1, 2);
    nc_def_var_deflate(file->id(), layout.U3_init->id(), 1, 1, 2);
    nc_def_var_deflate(file->id(), layout.U1->id(), 1, 1, 2);
    nc_def_var_deflate(file->id(), layout.U2->id(), 1, 1, 2);
    nc_def_var_deflate(file->id(), layout.U3->id(), 1, 1, 2);
#endif
    file->sync();
}

void KPNetCDFFile::readVars(HONEI_UNUSED KPInitialConditions& init) {
    //Create initial condidion variables
    layout.vars.init_B = file->get_var("init_B");

    //Create the timestep variables
    layout.vars.U1 = file->get_var("U1");
    layout.vars.U2 = file->get_var("U2");
    layout.vars.U3 = file->get_var("U3");
}


void KPNetCDFFile::writeAtts(HONEI_UNUSED const KPInitialConditions& init) {
    stringstream tmp;

    tmp.str("");
    time_t foo = time(NULL);
    tm* t = localtime(&foo);
    tmp << setfill('0');
    tmp << (1900+t->tm_year) << "-" << setw(2) << t->tm_mon << "-" << setw(2) << t->tm_mday;
    tmp << " " << setw(2) << t->tm_hour << ":" << setw(2) << t->tm_min << ":" << setw(2) << t->tm_sec;
    file->add_att("created_on", tmp.str().c_str());

    file->add_att("created_with", "coffee and salad");

    file->sync();
}


void KPNetCDFFile::readAtts(HONEI_UNUSED KPInitialConditions& init) {
    cout << "created_by   = " << file->get_att("created_by")->as_string(0) << endl;
    cout << "created_on   = " << file->get_att("created_on")->as_string(0) << endl;
    cout << "created_with = " << file->get_att("created_with")->as_string(0) << endl;
}


void KPNetCDFFile::writeInitialConditions(const KPInitialConditions& init) {
    //Write all attributes that describe the simulation
    file->add_att("B_location", init.getBLocation());
    file->add_att("nx", (int) init.getNx());
    file->add_att("ny", (int) init.getNy());
    file->add_att("dx", init.getDx());
    file->add_att("dy", init.getDy());
    file->add_att("g", init.getG());
    file->add_att("n", init.getN());
    file->add_att("dt_scale", init.getDtScale());
    file->add_att("desingularization_eps", init.getDesingularizationEps());
    file->add_att("time_integrator", init.getTimeIntegrator());
    file->add_att("bc_north", init.getBC().getNorth());
    file->add_att("bc_south", init.getBC().getSouth());
    file->add_att("bc_east", init.getBC().getEast());
    file->add_att("bc_west", init.getBC().getWest());
    file->add_att("bc_north_arg", init.getBC().getNorthArg());
    file->add_att("bc_south_arg", init.getBC().getSouthArg());
    file->add_att("bc_east_arg", init.getBC().getEastArg());
    file->add_att("bc_west_arg", init.getBC().getWestArg());

    //Write initial conditions
    layout.vars.init_B->put(init.getB(), init.getBNy(), init.getBNx());
    file->sync();
}


void KPNetCDFFile::readInitialConditions(KPInitialConditions& init,
        shared_ptr<Field>& init_B,
        shared_ptr<Field>& init_U1,
        shared_ptr<Field>& init_U2,
        shared_ptr<Field>& init_U3) {
    //Write all attributes that describe the simulation
    KPInitialConditions::DATA_LOCATION B_location = (KPInitialConditions::DATA_LOCATION) file->get_att("B_location")->as_int(0);
    unsigned int nx = file->get_att("nx")->as_int(0);
    unsigned int ny = file->get_att("ny")->as_int(0);
    float dx = file->get_att("dx")->as_float(0);
    float dy = file->get_att("dy")->as_float(0);
    //float g = file->get_att("g")->as_float(0);
    //float n = file->get_att("n")->as_float(0);
    //float dt_scale = file->get_att("dt_scale")->as_float(0);
    //float desingularization_eps = file->get_att("desingularization_eps")->as_float(0);
    KPInitialConditions::TIME_INTEGRATOR time_integrator = (KPInitialConditions::TIME_INTEGRATOR) file->get_att("time_integrator")->as_int(0);
    KPBoundaryConditions::TYPE bc_north = (KPBoundaryConditions::TYPE) file->get_att("bc_north")->as_int(0);
    KPBoundaryConditions::TYPE bc_south = (KPBoundaryConditions::TYPE) file->get_att("bc_south")->as_int(0);
    KPBoundaryConditions::TYPE bc_east = (KPBoundaryConditions::TYPE) file->get_att("bc_east")->as_int(0);
    KPBoundaryConditions::TYPE bc_west = (KPBoundaryConditions::TYPE) file->get_att("bc_west")->as_int(0);
    float bc_north_arg = file->get_att("bc_north_arg")->as_float(0);
    float bc_south_arg = file->get_att("bc_south_arg")->as_float(0);
    float bc_east_arg = file->get_att("bc_east_arg")->as_float(0);
    float bc_west_arg = file->get_att("bc_west_arg")->as_float(0);

    unsigned int BNx = nx + (B_location == KPInitialConditions::GRID_CELL_INTERSECTION);
    unsigned int BNy = ny + (B_location == KPInitialConditions::GRID_CELL_INTERSECTION);

    //Allocate for B, U[123]
    init_B.reset(new Field(BNx, BNy));
    init_U1.reset(new Field(nx, ny));
    init_U2.reset(new Field(nx, ny));
    init_U3.reset(new Field(nx, ny));

    //Get initial conditions
    layout.vars.init_B->get(init_B->data, init_B->ny, init_B->nx);
    layout.vars.U1->set_cur(nt);
    layout.vars.U1->get(init_U1->data, 1, init_U1->ny, init_U1->nx);
    layout.vars.U2->set_cur(nt);
    layout.vars.U2->get(init_U2->data, 1, init_U2->ny, init_U2->nx);
    layout.vars.U3->set_cur(nt);
    layout.vars.U3->get(init_U3->data, 1, init_U3->ny, init_U3->nx);

    //Create initial conditions variable
    KPBoundaryConditions bc(bc_north, bc_south, bc_east, bc_west,
            bc_north_arg, bc_south_arg, bc_east_arg, bc_west_arg);
    init = KPInitialConditions(nx, ny, dx, dy, bc, time_integrator,
            B_location, init_B->data,
            KPInitialConditions::GRID_CELL_CENTER, init_U1->data, init_U2->data, init_U3->data);
}


void KPNetCDFFile::writeTimeStep(shared_ptr<TimeStep> ts, int index) {
    t_index = index;
    if (t_index < 0) t_index = nt;

    float t = ts->time + time_offset;
    layout.vars.t->set_cur(t_index);
    layout.vars.t->put(&t, 1);

    layout.vars.U1->set_cur(t_index, 0, 0);
    layout.vars.U1->put(ts->U[0]->data, 1, ts->ny, ts->nx);

    layout.vars.U2->set_cur(t_index, 0, 0);
    layout.vars.U2->put(ts->U[1]->data, 1, ts->ny, ts->nx);

    layout.vars.U3->set_cur(t_index, 0, 0);
    layout.vars.U3->put(ts->U[2]->data, 1, ts->ny, ts->nx);

    file->sync();
    nt++;
}


void KPNetCDFFile::readTimeStep(shared_ptr<TimeStep> ts, int index) {
    if (index > 0) t_index = index;
    if (t_index >= nt) t_index = 0;

    layout.vars.t->set_cur(t_index);
    layout.vars.t->get(&ts->time, 1);

    layout.vars.U1->set_cur(t_index, 0, 0);
    layout.vars.U1->get(ts->U[0]->data, 1, ts->ny, ts->nx);

    layout.vars.U2->set_cur(t_index, 0, 0);
    layout.vars.U2->get(ts->U[1]->data, 1, ts->ny, ts->nx);

    layout.vars.U3->set_cur(t_index, 0, 0);
    layout.vars.U3->get(ts->U[2]->data, 1, ts->ny, ts->nx);

    file->sync();
    t_index++;
}
