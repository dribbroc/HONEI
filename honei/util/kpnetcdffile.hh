/* vim: set sw=4 sts=4 et foldmethod=syntax : */
/*
 * KPNetCDFFile.h
 *
 *  Created on: Apr 9, 2010
 *      Author: babrodtk
 */

#ifndef KPNETCDFFILE_H_
#define KPNETCDFFILE_H_

#include <honei/util/netcdf_datatypes.hh>
#include <honei/util/instantiation_policy.hh>

#include <tr1/memory>

#include <netcdfcpp.h>
#include <netcdf.h>

class KPNetCDFFile :
    public honei::InstantiationPolicy<KPNetCDFFile, honei::NonCopyable>
{
    public:
        KPNetCDFFile(std::string filename, const KPInitialConditions& init);
        KPNetCDFFile(std::string filename, KPInitialConditions& init,
                std::tr1::shared_ptr<Field>& init_B,
                std::tr1::shared_ptr<Field>& init_U1,
                std::tr1::shared_ptr<Field>& init_U2,
                std::tr1::shared_ptr<Field>& init_U3);
        ~KPNetCDFFile();

        void readInitialConditions(KPInitialConditions& init,
                std::tr1::shared_ptr<Field>& init_B,
                std::tr1::shared_ptr<Field>& init_U1,
                std::tr1::shared_ptr<Field>& init_U2,
                std::tr1::shared_ptr<Field>& init_U3);
        void writeInitialConditions(const KPInitialConditions& init);

        void writeTimeStep(std::tr1::shared_ptr<TimeStep> ts, int index=-1);
        void readTimeStep(std::tr1::shared_ptr<TimeStep> ts, int index=-1);

    private:
        void writeDims(const KPInitialConditions& init);
        void writeAtts(const KPInitialConditions& init);
        void writeVars(const KPInitialConditions& init);

        void readDims(KPInitialConditions& init);
        void readAtts(KPInitialConditions& init);
        void readVars(KPInitialConditions& init);

    private:
        std::tr1::shared_ptr<NcFile> file;

        struct {
            struct {
                NcDim* i;       //!< number of grid cell intersections
                NcDim* j;       //!< number of grid cell intersections
                NcDim* x;       //!< number of grid cells
                NcDim* y;       //!< number of grid cells
                NcDim* t;       //!< Time
            } dims;

            struct {
                NcVar* init_B;   //!< Initial bathymetry
                NcVar* U1;       //!< U1 = w
                NcVar* U2;       //!< U2 = hu
                NcVar* U3;       //!< U3 = hv
                NcVar* i;        //!< i
                NcVar* j;        //!< j
                NcVar* x;        //!< x
                NcVar* y;        //!< y
                NcVar* t;        //!< time
            } vars;

            struct {
                NcAtt* i;
                NcAtt* j;
                NcAtt* x;
                NcAtt* y;
                NcAtt* t;
            } atts;
        } layout;
        long nt;
        long t_index;
        float time_offset;
};

#endif /* KPNETCDFFILE_H_ */
