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
 * \operand i Alignment offset to the first element of the first block of first matrix.
 * \operand j Number of columns of second matrix.
 * \operand k Number of DMA Lists for the result matrix.
 * \operand l Pointer to DMA transfer list pointers for the result matrix.
 * \operand m Pointer to DMA transfer list sizes for the result matrix.
 * \operand n Pointer to DMA transfer list effective addresses for the result matrix.
 * \operand o Number of columns of first matrix.
 */
void product_dense_matrix_dense_matrix_float(const Instruction & inst)
{
    debug_value(0);
    EffectiveAddress ea_a(inst.a.ea);

    Allocation * block_a[2] = { acquire_block(), acquire_block() };
    Allocation * block_b[2] = { acquire_block(), acquire_block() };
    Allocation * block_r[2] = { acquire_block(), acquire_block() };

    unsigned nr_of_lists(inst.g.u); // Number of Lists to be double buffered
    unsigned a_cols(inst.o.u);
    unsigned a_t_size(inst.size);
    unsigned a_typed_offset(inst.i.u); // if we are in lower part of A our first row may have an offset
    unsigned long long list_sizes[nr_of_lists] __attribute__((aligned(16)));
    unsigned long long list_eahs[nr_of_lists] __attribute__((aligned(16)));
    EffectiveAddress list_ptrs[nr_of_lists] __attribute__((aligned(16)));

    mfc_get(&list_sizes, inst.e.ea, multiple_of_sixteen(nr_of_lists * 8), 3, 0, 0);
    mfc_get(&list_ptrs, inst.d.ea, multiple_of_sixteen(nr_of_lists * 8), 3, 0, 0);
    mfc_get(&list_eahs, inst.f.ea, multiple_of_sixteen(nr_of_lists * 8), 3, 0, 0);
    mfc_write_tag_mask(1 << 3);
    mfc_read_tag_status_all();
    debug_value(1);
    unsigned a_counter(inst.b.u);
    unsigned a_size(a_counter > 1 ? a_t_size : multiple_of_sixteen(inst.c.u));
    unsigned a_current(0), a_next(1);
    unsigned b_current(0), b_next(1);
    unsigned r_current(0), r_next(1);

    Allocation * block_list[2] = { acquire_block(), acquire_block() };
    Pointer<ListElement> list_ptr[2] = { { block_list[0]->address }, { block_list[1]->address } };

    unsigned b_counter(0);
    debug_get(list_ptrs[0], list_ptr[b_current].untyped, multiple_of_sixteen(sizeof(ListElement) * list_sizes[0]));
    mfc_get(list_ptr[b_current].untyped, list_ptrs[0], multiple_of_sixteen(sizeof(ListElement) * list_sizes[0]), 1, 0, 0);
    mfc_write_tag_mask(1 << 1);
    mfc_read_tag_status_all();
    debug_value(2);
    Pointer<float> a[2] = { { block_a[0]->address }, { block_a[1]->address } };
    Pointer<float> b[2] = { { block_b[0]->address }, { block_b[1]->address } };
    Pointer<float> r[2] = { { block_r[0]->address }, { block_r[1]->address } };
    fill(r[0].untyped, 16384, 0.0f);
    unsigned b_size(list_sizes[b_counter]);
    debug_getl(list_eahs[b_counter], b[b_current].untyped, list_sizes[0] * sizeof(ListElement));
    mfc_getl(b[b_current].untyped, list_eahs[b_counter], list_ptr[0].untyped, list_sizes[0] * sizeof(ListElement), 6, 0, 0);
    debug_value(3);

    unsigned a_nextsize, b_nextsize;

    debug_get(ea_a, a[a_current].untyped, a_size);
    debug_value(33);
    mfc_get(a[a_current].untyped, ea_a, a_size, a_current, 0, 0);

    ea_a += a_size;
    debug_value(4);

    unsigned r_nr_of_lists(inst.k.u); // Number of Lists to be putted
    unsigned long long r_list_sizes[r_nr_of_lists] __attribute__((aligned(16)));
    unsigned long long r_list_eahs[r_nr_of_lists] __attribute__((aligned(16)));
    EffectiveAddress r_list_ptrs[r_nr_of_lists] __attribute__((aligned(16)));
    debug_value(8);
    mfc_get(&r_list_sizes, inst.m.ea, multiple_of_sixteen(r_nr_of_lists * 8), 5, 0, 0);
    mfc_get(&r_list_ptrs, inst.l.ea, multiple_of_sixteen(r_nr_of_lists * 8), 5, 0, 0);
    mfc_get(&r_list_eahs, inst.n.ea, multiple_of_sixteen(r_nr_of_lists * 8), 5, 0, 0);
    mfc_write_tag_mask(1 << 5);
    mfc_read_tag_status_all();
    debug_value(9);
    Allocation * r_block_list[2] = { acquire_block(), acquire_block() };
    Pointer<ListElement> r_list_ptr[2] = { { r_block_list[0]->address }, { r_block_list[1]->address } };
    unsigned r_elem_t_size = ((inst.j.u * 4) + 16);
    union lsaddr
    {
        void * ptr;
        unsigned long long value;
    };
    lsaddr lsa = { r[r_current].untyped };
    unsigned act_list(0);

    unsigned long b_offset(inst.h.u / 4); // Start offset of the row with index 1 (!) (row 0 has always index 0)
    unsigned b_vecs(inst.j.u / 4);
    unsigned ar(0);

    while (a_counter > 1) // db for A
    {
        debug_value(22);
        a_nextsize = (a_counter == 2 ? multiple_of_sixteen(inst.c.u) : a_t_size);

        debug_get(ea_a, a[a_next].untyped, a_nextsize);
        mfc_get(a[a_next].untyped, ea_a, a_nextsize, a_next, 0, 0);

        ea_a += a_nextsize;

        unsigned get_counter(0);
        unsigned long a_elem(0); // The actual considered element of matrix a
        unsigned a_rows(a_size / 4 / a_cols);
        unsigned r_offset(0);
        fill(r[r_current].untyped, 16384, 0.0f);

        for( ; ar < a_rows ; ar++)
        {
            b_counter = get_counter;
            unsigned e_offset(0);

            while (b_counter <= nr_of_lists - 1) // db for Lists of B
            {
                debug_value(23);
                b_counter++;
                get_counter = b_counter % nr_of_lists;
                b_nextsize = (list_sizes[get_counter]);

                debug_get(list_ptrs[get_counter], list_ptr[b_next].untyped, multiple_of_sixteen(sizeof(ListElement) * list_sizes[get_counter]));
                mfc_get(list_ptr[b_next].untyped, list_ptrs[get_counter], multiple_of_sixteen(sizeof(ListElement) * list_sizes[get_counter]), 4, 0, 0);
                mfc_write_tag_mask(1 << 4);
                mfc_read_tag_status_all();

                debug_getl(list_eahs[get_counter], b[b_next].untyped, list_sizes[get_counter] * sizeof(ListElement));
                mfc_getl(b[b_next].untyped, list_eahs[get_counter], list_ptr[b_next].untyped, list_sizes[get_counter] * sizeof(ListElement), 9 + b_next, 0, 0);

                mfc_write_tag_mask(1 << 9 + b_current);
                mfc_read_tag_status_all();

                unsigned b_nr_rows(list_sizes[b_counter-1]);
                unsigned b_vec_idx(0);

                for (unsigned br(0) ; br < b_nr_rows ; br++, a_elem++)
                {
                    unsigned r_idx(ar * (b_vecs + 1));

                    for(unsigned i(0) ; i < b_vecs ; i++)
                    {
                        vector float temp = b[b_current].vectorised[b_vec_idx]; // temp version needed, cause original matrix must not be changed!
                        extract(temp, b[b_current].vectorised[b_vec_idx + 1], e_offset);

                        vector float r_temp(r[r_current].vectorised[r_idx]);
                        extract(r_temp, r[r_current].vectorised[r_idx+1], r_offset);
                        r_temp = spu_madd(spu_splats(a[a_current].typed[a_elem + a_typed_offset]), temp, r_temp);
                        insert(r[r_current].vectorised[r_idx], r[r_current].vectorised[r_idx + 1], r_temp, r_offset);

                        //Subscriptable<float> fv = { temp };
                        //Subscriptable<float> rs = { r[r_current].vectorised[r_idx] };
                        //printf("Multipyling: a_elem: %f * %f %f %f %f \n", a[a_current].typed[a_elem+ a_typed_offset], fv.array[0], fv.array[1], fv.array[2], fv.array[3]);
                        //printf("Result an pos: %u:     %f %f %f %f \n", r_idx, rs.array[0], rs.array[1], rs.array[2], rs.array[3]);

                        b_vec_idx++;
                        r_idx++;
                        if (i == b_vecs - 1 && e_offset > 0)
                        {
                            b_vec_idx++;
                        }
                    }

                    e_offset = (e_offset + b_offset) % 4;
                }

                unsigned b_temp(b_next);
                b_next = b_current;
                b_current = b_temp;
                b_size = b_nextsize;

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
            act_list++;
            lsa.value += r_list_sizes[act_list] * r_elem_t_size;
        }

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
    unsigned get_counter(0);
    unsigned long a_elem(0); // The actual considered element of matrix a
    unsigned a_rows(a_size / 4 / a_cols);
    unsigned r_offset(0);
    fill(r[r_current].untyped, 16384, 0.0f);

    for( ; ar < a_rows ; ar++)
    {
        b_counter = get_counter;
        unsigned e_offset(0);

        while (b_counter <= nr_of_lists - 1) // db for Lists of B
        {
            debug_value(55);
            b_counter++;
            get_counter = b_counter % nr_of_lists;
            debug_value(get_counter);
            b_nextsize = (list_sizes[get_counter]);

            debug_get(list_ptrs[get_counter], list_ptr[b_next].untyped, multiple_of_sixteen(sizeof(ListElement) * list_sizes[get_counter]));
            mfc_get(list_ptr[b_next].untyped, list_ptrs[get_counter], multiple_of_sixteen(sizeof(ListElement) * list_sizes[get_counter]), 4, 0, 0);
            mfc_write_tag_mask(1 << 4);
            mfc_read_tag_status_all();

            debug_getl(list_eahs[get_counter], b[b_next].untyped, list_sizes[get_counter] * sizeof(ListElement));
            mfc_getl(b[b_next].untyped, list_eahs[get_counter], list_ptr[b_next].untyped, list_sizes[get_counter] * sizeof(ListElement), 9 + b_next, 0, 0);

            mfc_write_tag_mask(1 << 9 + b_current);
            mfc_read_tag_status_all();

            unsigned b_nr_rows(list_sizes[b_counter-1]);
            unsigned b_vec_idx(0);

            for (unsigned br(0) ; br < b_nr_rows ; br++, a_elem++)
            {
                unsigned r_idx(ar * (b_vecs + 1));

                for(unsigned i(0) ; i < b_vecs ; i++)
                {
                    vector float temp = b[b_current].vectorised[b_vec_idx]; // temp version needed, cause original matrix must not be changed!
                    extract(temp, b[b_current].vectorised[b_vec_idx + 1], e_offset);

                    vector float r_temp(r[r_current].vectorised[r_idx]);
                    extract(r_temp, r[r_current].vectorised[r_idx+1], r_offset);
                    r_temp = spu_madd(spu_splats(a[a_current].typed[a_elem + a_typed_offset]), temp, r_temp);
                    insert(r[r_current].vectorised[r_idx], r[r_current].vectorised[r_idx + 1], r_temp, r_offset);
/*
                    Subscriptable<float> fv = { temp };
                    Subscriptable<float> rs = { r_temp };
                    printf("Multipyling: a_elem: %f * %f %f %f %f \n", a[a_current].typed[a_elem+ a_typed_offset], fv.array[0], fv.array[1], fv.array[2], fv.array[3]);
                    printf("Result an pos: %u:     %f %f %f %f \n", r_idx, rs.array[0], rs.array[1], rs.array[2], rs.array[3]);
*/
                    b_vec_idx++;
                    r_idx++;
                    if (i == b_vecs - 1 && e_offset > 0)
                    {
                        b_vec_idx++;
                    }
                }

                e_offset = (e_offset + b_offset) % 4;
            }

            unsigned b_temp(b_next);
            b_next = b_current;
            b_current = b_temp;
            b_size = b_nextsize;

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
            mfc_putl(lsa.ptr, r_list_eahs[act_list], r_list_ptr[r_current].untyped, r_list_sizes[act_list] * sizeof(ListElement), 6 + r_next, 0, 0);
            act_list++;
            lsa.value += r_list_sizes[act_list] * r_elem_t_size;
        }

    mfc_write_tag_mask(1 << 6);
    mfc_read_tag_status_all();

    release_all_blocks();

}
