program feat
real*8 , dimension(0:4) :: a
real*8 , dimension(0:4) :: b
integer , dimension(0:7) :: border_mask
integer*8  size;
real*8 tolerance;
a(0) = 1;
a(1) = 2;
a(2) = 3;
a(3) = 4;
a(4) = 5;
b(0) = 2;
b(1) = 4;
b(2) = 6;
b(3) = 8;
b(4) = 10;
size = 5;
print *, 'a=', a;
print *, 'b=',  b;
call sum(a, b, size);
print *, 'result=', a;

!!!!!!!!!!!!!!!!!!!!! MG krams
call is_smoother(.FALSE.);
border_mask(0) = 2;
border_mask(1) = 2;
border_mask(2) = 2;
border_mask(3) = 2;
border_mask(4) = 2;
border_mask(5) = 1;
border_mask(6) = 2;
border_mask(7) = 2;
call set_border_mask(border_mask);
call set_lvl_range(1, 7);
call set_n_max_iter(16);
call set_initial_zero(.FALSE.);
tolerance = 1e-8;
call set_tolerance(tolerance);
call set_convergence_check(.TRUE.);
call set_n_pre_smooth(2);
call set_n_post_smooth(2);
!!! gelitten! rechne doch selber aus
call set_n_max_iter_coarse(4711);
tolerance = 1e-2;
call set_tolerance_coarse(tolerance);
tolerance = 1.
call set_adapt_correction_factor(tolerance);
!!call push_rhs(a, size);


call run();
end program feat
