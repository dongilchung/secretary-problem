function u_x = cal_u_x( power_params, reference, x_vector )


try
[r,c] = size( x_vector );

if r > c
    x_vector = x_vector';
end

if length( power_params ) == 1
    c = power_params;
    b = power_params;
    a = 1;
elseif length( power_params ) == 2
    c = power_params(1);
    b = power_params(1);
    a = power_params(2);
else
    c = power_params(3);
    b = power_params(2);
    a = power_params(1);
end

x_vector = [0, x_vector, 150];

diff_x = x_vector - reference;
abs_diff_x = abs(diff_x);
abs_power_x( x_vector < reference ) = abs_diff_x( x_vector < reference ).^b;
abs_power_x( x_vector >= reference ) = abs_diff_x( x_vector >= reference ).^c;
abs_power_x( x_vector < reference ) = -a*abs_power_x( x_vector < reference );

u_x = (abs_power_x - abs_power_x(1))/(abs_power_x(end) - abs_power_x(1));
u_x = u_x(2:end-1);

catch me

    x_vector
    abs_power_x
    reference

end