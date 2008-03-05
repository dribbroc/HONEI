/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Sven Mallach <sven.mallach@honei.org>
 *
 * This file is part of the LA C++ library. LibLa is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibLa is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <honei/cell/cell.hh>
#include <honei/cell/libutil/allocator.hh>
#include <honei/cell/libutil/transfer.hh>

#include <spu_intrinsics.h>
#include <spu_mfcio.h>

using namespace honei::cell;

/*
 * dense_dense_float_matrix_product
 *
 * Calculate the product of two dense matrices.
 *
 * \size default transfer size for first matrix.
 * \operand a Base address of first matrix.
 * \operand b Number of transfers for first matrix.
 * \operand c Size of the last transfers for the first matrix (use multiple_of_sixteen())
 * \operand d Pointer to DMA transfer list pointers for the second matrix.
 * \operand e Pointer to DMA transfer list sizes for the second matrix.
 * \operand f Pointer to DMA transfer list effective addresses for the second matrix.
 * \operand g Number of DMA Lists for the second matrix.
 * \operand h Alignment offset of the second row of the second matrix.
 * \operand i Number of columns of first matrix.
 * \operand j Number of columns of second matrix.
 * \operand k Number of DMA Lists for the result matrix.
 * \operand l Pointer to DMA transfer list pointers for the result matrix.
 * \operand m Pointer to DMA transfer list sizes for the result matrix.
 * \operand n Pointer to DMA transfer list effective addresses for the result matrix.
 */
void product_dense_matrix_dense_matrix_float(const Instruction & inst)
{
    debug_value(0);

    const unsigned nr_of_lists(inst.g.u); // Number of Lists to be double buffered
    const unsigned a_cols(inst.i.u);
    const unsigned a_t_size(inst.size);
    unsigned list_sizes[nr_of_lists] __attribute__((aligned(16)));
    unsigned long long list_eahs[nr_of_lists] __attribute__((aligned(16)));
    EffectiveAddress list_ptrs[nr_of_lists] __attribute__((aligned(16)));

    // GET list information for Matrix B
    mfc_get(&list_sizes, inst.e.ea, multiple_of_sixteen(nr_of_lists * 4), 3, 0, 0);
    mfc_get(&list_ptrs, inst.d.ea, multiple_of_sixteen(nr_of_lists * 8), 3, 0, 0);
    mfc_get(&list_eahs, inst.f.ea, multiple_of_sixteen(nr_of_lists * 8), 3, 0, 0);

    debug_value(1);
    unsigned a_counter(inst.b.u), b_counter(0);
    unsigned a_size(a_counter > 1 ? a_t_size : multiple_of_sixteen(inst.c.u)), a_nextsize;
    unsigned a_current(0), a_next(1);
    unsigned b_current(0), b_next(1);
    unsigned r_current(0), r_next(1);

    Allocation * block_list[2] = { acquire_block(), acquire_block() };
    Pointer<ListElement> list_ptr[2] = { { block_list[0]->address }, { block_list[1]->address } };

    // Assure that list information has arrived.
    mfc_write_tag_mask(1 << 3);
    mfc_read_tag_status_all();

    // GET the first DMA List for Matrix B
    debug_get(list_ptrs[0], list_ptr[b_current].untyped, multiple_of_sixteen(sizeof(ListElement) * list_sizes[0]));
    mfc_get(list_ptr[b_current].untyped, list_ptrs[0], multiple_of_sixteen(sizeof(ListElement) * list_sizes[0]), 10, 0, 0);
    debug_value(2);

    const unsigned r_nr_of_lists(inst.k.u); // Number of Lists to be putted
    unsigned r_list_sizes[r_nr_of_lists] __attribute__((aligned(16)));
    unsigned long long r_list_eahs[r_nr_of_lists] __attribute__((aligned(16)));
    EffectiveAddress r_list_ptrs[r_nr_of_lists] __attribute__((aligned(16)));

    Allocation * block_a[2] = { acquire_block(), acquire_block() };
    Allocation * block_b[4] = { acquire_block(), acquire_block(), acquire_block(), acquire_block() };
    Allocation * block_r[2] = { acquire_block(), acquire_block() };

    Pointer<float> a[2] = { { block_a[0]->address }, { block_a[1]->address } };
    Pointer<float> b[2] = { { block_b[0]->address }, { block_b[2]->address } };
    Pointer<float> r[2] = { { block_r[0]->address }, { block_r[1]->address } };

    // Assure that first DMA List (B) has arrived
    mfc_write_tag_mask(1 << 10);
    mfc_read_tag_status_all();

    // GET(L) the first DMA transfers using the first List for B
    debug_getl(list_eahs[b_counter], b[b_current].untyped, list_sizes[0] * sizeof(ListElement));
    mfc_getl(b[b_current].untyped, list_eahs[b_counter], list_ptr[0].untyped, list_sizes[0] * sizeof(ListElement), 6, 0, 0);
    debug_value(3);

    // GET the next DMA list for B
    debug_get(list_ptrs[1 % nr_of_lists], list_ptr[1].untyped, multiple_of_sixteen(sizeof(ListElement) * list_sizes[1 % nr_of_lists]));
    mfc_get(list_ptr[1].untyped, list_ptrs[1 % nr_of_lists], multiple_of_sixteen(sizeof(ListElement) * list_sizes[1 % nr_of_lists]), 10, 0, 0);

    fill(r[r_current].untyped, 16384, 0.0f);

    EffectiveAddress ea_a(inst.a.ea);
    // GET the first part of Matrix A
    debug_get(ea_a, a[a_current].untyped, a_size);
    mfc_get(a[a_current].untyped, ea_a, a_size, a_current, 0, 0);
    ea_a += a_size;

    // GET list information for Matrix R
    mfc_get(&r_list_sizes, inst.m.ea, multiple_of_sixteen(r_nr_of_lists * 4), 5, 0, 0);
    mfc_get(&r_list_ptrs, inst.l.ea, multiple_of_sixteen(r_nr_of_lists * 8), 5, 0, 0);
    mfc_get(&r_list_eahs, inst.n.ea, multiple_of_sixteen(r_nr_of_lists * 8), 5, 0, 0);

    debug_value(4);

    Allocation * r_block_list[2] = { acquire_block(), acquire_block() };
    Pointer<ListElement> r_list_ptr[2] = { { r_block_list[0]->address }, { r_block_list[1]->address } };
    const unsigned r_elem_t_size((inst.j.u * 4) + 16);
    union lsaddr
    {
        void * ptr;
        unsigned long long value;
    };
    lsaddr lsa = { r[r_current].untyped };
    unsigned act_list(0);

    const unsigned long b_offset(inst.h.u / 4); // Start offset of the row with index 1 (!) (row 0 has always index 0)
    const unsigned b_vecs(inst.j.u / 4);
    unsigned ar(0);
    unsigned r_offset(0);

    // Assure that list information for Matrix R has arrived
    mfc_write_tag_mask(1 << 5);
    mfc_read_tag_status_all();

    debug_value(5);

    while (a_counter > 1) // db for A
    {
        debug_value(22);
        a_nextsize = (a_counter == 2 ? multiple_of_sixteen(inst.c.u) : a_t_size);

        debug_get(ea_a, a[a_next].untyped, a_nextsize);
        mfc_get(a[a_next].untyped, ea_a, a_nextsize, a_next, 0, 0);

        ea_a += a_nextsize;

        unsigned get_counter(0), get_next_counter(0);
        unsigned long a_elem(0); // The actual considered element of matrix a
        unsigned a_rows(a_size / 4 / a_cols);

        mfc_write_tag_mask(1 << a_current);
        mfc_write_tag_mask(1 << 6 + r_current);
        mfc_read_tag_status_all();

        for( ; ar < a_rows ; ar++)
        {
            b_counter = get_counter;
            unsigned e_offset(0);

            while (b_counter <= nr_of_lists - 1) // db for Lists of B
            {
                debug_value(23);
                b_counter++;
                get_counter = b_counter % nr_of_lists;
                get_next_counter = (b_counter + 1) % nr_of_lists;

                // Assure that next DMA List is ready
                mfc_write_tag_mask(1 << 10);
                mfc_read_tag_status_all();

                // GETL next
                debug_getl(list_eahs[get_counter], b[b_next].untyped, list_sizes[get_counter] * sizeof(ListElement));
                mfc_getl(b[b_next].untyped, list_eahs[get_counter], list_ptr[b_next].untyped, list_sizes[get_counter] * sizeof(ListElement), 9 + b_next, 0, 0);

                // Assure that last DMA GETL has finished
                mfc_write_tag_mask(1 << 9 + b_current);
                mfc_read_tag_status_all();

                // GET second next DMA List for the following iteration
                debug_get(list_ptrs[get_next_counter], list_ptr[b_current].untyped, multiple_of_sixteen(sizeof(ListElement) * list_sizes[get_next_counter]));
                mfc_get(list_ptr[b_current].untyped, list_ptrs[get_next_counter], multiple_of_sixteen(sizeof(ListElement) * list_sizes[get_next_counter]), 10, 0, 0);

                const unsigned b_nr_rows(list_sizes[b_counter - 1]);
                unsigned b_vec_idx(0);

                vector float r_temps[b_vecs];
                const unsigned r_idx(ar * (b_vecs + 1));

                for (unsigned i(0) ; i < b_vecs ; i++)
                {
                    r_temps[i] = r[r_current].vectorised[r_idx + i];
                    extract(r_temps[i], r[r_current].vectorised[r_idx + i + 1], r_offset);
                }

                for (unsigned br(0) ; br < b_nr_rows ; br++, a_elem++)
                {
                    for(unsigned i(0) ; i < b_vecs ; i++, b_vec_idx++)
                    {
                        vector float temp = b[b_current].vectorised[b_vec_idx]; // temp version needed, cause original matrix must not be changed!
                        extract(temp, b[b_current].vectorised[b_vec_idx + 1], e_offset);

                        r_temps[i] = spu_madd(spu_splats(a[a_current].typed[a_elem]), temp, r_temps[i]);
                    }

                    b_vec_idx++;
                    e_offset = (e_offset + b_offset) % 4;
                }

                for (unsigned i(0) ; i < b_vecs ; i++)
                {
                    insert(r[r_current].vectorised[r_idx + i], r[r_current].vectorised[r_idx + i + 1], r_temps[i], r_offset);
                }

                unsigned b_temp(b_next);
                b_next = b_current;
                b_current = b_temp;

            } // end while B

            r_offset = (r_offset + b_offset) % 4;
        } // end for

        // PUT a_rows
        unsigned iterations(0);
        unsigned size_acc(0);
        while (size_acc != a_rows)
        {
            size_acc += r_list_sizes[act_list + iterations];
            iterations++;
        }
        for (unsigned x(0) ; x < iterations ; x++)
        {
            debug_value(1441);
            debug_value(r_current);
            debug_value(act_list);

            debug_get(r_list_ptrs[act_list], r_list_ptr[r_current].untyped, multiple_of_sixteen(sizeof(ListElement) * r_list_sizes[act_list]));
            mfc_get(r_list_ptr[r_current].untyped, r_list_ptrs[act_list], multiple_of_sixteen(sizeof(ListElement) * r_list_sizes[act_list]), 4, 0, 0);

            mfc_write_tag_mask(1 << 4);
            mfc_read_tag_status_all();

            debug_putl(r_list_eahs[act_list], lsa.ptr, r_list_sizes[act_list] * sizeof(ListElement));
            mfc_putl(lsa.ptr, r_list_eahs[act_list], r_list_ptr[r_current].untyped, r_list_sizes[act_list] * sizeof(ListElement), 6 + r_next, 0, 0);

            lsa.value += r_list_sizes[act_list] * r_elem_t_size;
            act_list++;
        }

        fill(r[r_next].untyped, 16384, 0.0f);

        unsigned a_temp(a_next);
        a_next = a_current;
        a_current = a_temp;
        a_size = a_nextsize;
        a_counter--;
        ar = 0;
        unsigned r_temp(r_next);
        r_next = r_current;
        r_current = r_temp;
        lsa.ptr = r[r_current].untyped;

    } // end while A
    debug_value(44);
    mfc_write_tag_mask(1 << a_current);
    mfc_read_tag_status_all();
    unsigned get_counter(0), get_next_counter(0);
    unsigned long a_elem(0); // The actual considered element of matrix a
    unsigned a_rows(a_size / 4 / a_cols);

    for( ; ar < a_rows ; ar++)
    {
        b_counter = get_counter;
        unsigned e_offset(0);

        while (b_counter <= nr_of_lists - 1) // db for Lists of B
        {
            debug_value(55);
            b_counter++;
            get_counter = b_counter % nr_of_lists;
            get_next_counter = (b_counter + 1) % nr_of_lists;
            debug_value(get_counter);

            // Assure that next DMA List is ready
            mfc_write_tag_mask(1 << 10);
            mfc_read_tag_status_all();

            // GETL next
            debug_getl(list_eahs[get_counter], b[b_next].untyped, list_sizes[get_counter] * sizeof(ListElement));
            mfc_getl(b[b_next].untyped, list_eahs[get_counter], list_ptr[b_next].untyped, list_sizes[get_counter] * sizeof(ListElement), 9 + b_next, 0, 0);

            // Assure that last DMA GETL has finished
            mfc_write_tag_mask(1 << 9 + b_current);
            mfc_read_tag_status_all();

            // GET second next DMA List for the following iteration
            debug_get(list_ptrs[get_next_counter], list_ptr[b_current].untyped, multiple_of_sixteen(sizeof(ListElement) * list_sizes[get_next_counter]));
            mfc_get(list_ptr[b_current].untyped, list_ptrs[get_next_counter], multiple_of_sixteen(sizeof(ListElement) * list_sizes[get_next_counter]), 10, 0, 0);


            const unsigned b_nr_rows(list_sizes[b_counter - 1]);
            unsigned b_vec_idx(0);

            vector float r_temps[b_vecs];
            const unsigned r_idx(ar * (b_vecs + 1));

            for (unsigned i(0) ; i < b_vecs ; i++)
            {
                r_temps[i] = r[r_current].vectorised[r_idx + i];
                extract(r_temps[i], r[r_current].vectorised[r_idx + i + 1], r_offset);
            }

            for (unsigned br(0) ; br < b_nr_rows ; br++, a_elem++)
            {
                for(unsigned i(0) ; i < b_vecs ; i++, b_vec_idx++)
                {
                    vector float temp = b[b_current].vectorised[b_vec_idx]; // temp version needed, cause original matrix must not be changed!
                    extract(temp, b[b_current].vectorised[b_vec_idx + 1], e_offset);

                    r_temps[i] = spu_madd(spu_splats(a[a_current].typed[a_elem]), temp, r_temps[i]);
                }

                b_vec_idx++;
                e_offset = (e_offset + b_offset) % 4;
            }

            for (unsigned i(0) ; i < b_vecs ; i++)
            {
                insert(r[r_current].vectorised[r_idx + i], r[r_current].vectorised[r_idx + i + 1], r_temps[i], r_offset);
            }

            unsigned b_temp(b_next);
            b_next = b_current;
            b_current = b_temp;
        }

        r_offset = (r_offset + b_offset) % 4;

    }

    debug_value(7);

    // PUT a_rows
    unsigned iterations(0);
    unsigned size_acc(0);
    while (size_acc != a_rows)
    {
        size_acc += r_list_sizes[act_list + iterations];
        iterations++;
    }

    for (unsigned x(0) ; x < iterations ; x++)
    {
        debug_value(2882);
        debug_value(r_current);
        debug_value(act_list);
        debug_get(r_list_ptrs[act_list], r_list_ptr[r_current].untyped, multiple_of_sixteen(sizeof(ListElement) * r_list_sizes[act_list]));
        mfc_get(r_list_ptr[r_current].untyped, r_list_ptrs[act_list], multiple_of_sixteen(sizeof(ListElement) * r_list_sizes[act_list]), 4, 0, 0);
        mfc_write_tag_mask(1 << 4);
        mfc_read_tag_status_all();

        debug_putl(r_list_eahs[act_list], lsa.ptr, r_list_sizes[act_list] * sizeof(ListElement));
        mfc_putl(lsa.ptr, r_list_eahs[act_list], r_list_ptr[r_current].untyped, r_list_sizes[act_list] * sizeof(ListElement), 6, 0, 0);
        lsa.value += r_list_sizes[act_list] * r_elem_t_size;
        act_list++;
    }

    mfc_write_tag_mask(1 << 6);
    mfc_read_tag_status_all();

    release_all_blocks();

}
