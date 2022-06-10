---- Internal Gear Pump
-- 3D_it group
-- Mohammadmahdi Ataei
-- Arash Kadkhodaei Elyaderani
-- Principal Supervisor: Prof. Dr.-Ing. Stefan Scherbarth
-- Case Study Cyber Physical Production Systems using Additive manufacturing
-- Deggendorf Institute of Technology
-- Last edit: 12.06.2022


---- Input parameters using Interface
z_2 = ui_numberBox("Number Of Teeth for Rotor Gear", 25);  --number of teeth Rotor
z_1 = ui_numberBox("Number Of Teeth for Idler Gear", 17);  --number of teeth Idler
m = ui_numberBox("Module Of Gear", 4);  -- Module
alpha_t = ui_scalarBox("Pressure Angle", 20, 1);  -- Pressure angle
h_a_coef_p = ui_scalarBox("Addendum Coef(mm)", 1, 0.05);  -- Addendum height
h_f_coef_p = ui_scalarBox("Dedendum Coef(mm)", 1.25, 0.05);  -- Dedendum height
x_coef_int = ui_scalarBox("Internal Profile Shift(mm)", 0.1, 0.05);  -- Profile shift factor for internal gear
x_coef_ext = ui_scalarBox("External Profile Shift(mm)", 0, 0.05);  -- Profile shift factor for external gear
b = ui_numberBox("Width(mm)", 10);  -- Thickness of the gear
rotation = ui_numberBox("Rotate", 0);  -- Rotation


---- Function for involute profile points
function involute_t(r_b, inv_alpha) -- r_b
    return v(r_b * (math.sin(inv_alpha) - inv_alpha * math.cos(inv_alpha)), r_b * (math.cos(inv_alpha) + inv_alpha * math.sin(inv_alpha)))
end

---- Function for involute angle
function involute_angle(r_p1, r_p2)
    return (math.sqrt((r_p2 * r_p2 - r_p1 * r_p1) / (r_p1 * r_p1)))
end

---- Function for rotation matrix
function p_rotate(a, points) -- a = angle
    return v(math.cos(a) * points.x + math.sin(a) * points.y, math.cos(a) * points.y - math.sin(a) * points.x)
end

---- Function for mirroring
function mirror(point)
    return v(-point.x, point.y)
end

---- Function for  Internal and External Gear Profiles
function gearProfile(z, m_n, alpha_t, x_coef, h_a_coef, h_f_coef, b)
    local inv_xy = {}

    local alpha = alpha_t * math.pi / 180;  -- Pressure angle
    c = 0.1 * m_n;  -- Clearance

    -- x_coef is Profile shift
    if x_coef < 0 then
        x_coef = 0;
    end

    -- Formula calculations for the gear profile

    local d_p = z * m_n -- Pitch Diamter
    local r_p = d_p / 2 -- Pitch radius

    local d_b = d_p * math.cos(alpha);  -- Base diameter of gear
    r_b = d_b / 2 -- Base radius

    -- h_f_coef is dedendum coefficient
    local h_f = m_n * (h_f_coef) -- Dedendum

    -- h_a_coef is addendum coefficient
    local h_a = m_n * (h_a_coef);  -- Addendum
    local d_a = (d_p + (2 * m_n * (1 + x_coef)));
    r_a = d_a / 2;  -- Addendum radius

    local d_f = m_n * z + 2 * x_coef * m_n - 2 * (h_f);
    r_f = d_f / 2;  -- Root_radius

    -- Involute function
    inv_a = math.tan(alpha) - alpha;

    -- True Involute Diameter:
    local fp = ((m_n * ((math.pi / 4) - (math.tan(alpha))) - c * math.tan(alpha)) * (1 + math.sin(alpha))) / (math.cos(alpha))
    local d_TIF = math.sqrt(math.pow(d_p * math.sin(alpha) - 2 * (h_a - (m_n * x_coef) - fp * (1 - math.sin(alpha))) / (math.sin(alpha)), 2) + d_b * d_b) -- True involute diameter
    r_TIF = d_TIF / 2;  --radius

    -- Tooth thickness -->> pitch circle
    local S_0 = m_n * ((math.pi / 2) + 2 * x_coef * math.tan(alpha));

    -- Tooth thickness on the form circle
    local alpha_f     = math.acos((d_p * math.cos(alpha)) / d_TIF);  -- Involute angle along form circle
    local invFunction = math.tan(alpha_f) - alpha_f;  -- Involute function
    local s_f         = d_TIF * ((S_0 / d_p) + inv_a - invFunction) + (x_coef * math.tan(alpha_f));  -- Thickness of the gear teeth along form circle
    local angel       = (s_f) / (0.5 * d_TIF);  -- Angle swept along form circle for corresponding thickeness


    -- loops for generation points
    local steps = 30;
    for i = 1, z do
        local th = 2 * math.pi
        local AngelE = involute_angle(r_TIF, r_a);  -- End asngle
        for j = 1, steps do
            inv_xy[#inv_xy + 1] = p_rotate(th * i / z, involute_t(r_TIF, ((AngelE) * j / steps))) -- firsr involute face
        end
        for j = steps, 1, -1 do
            inv_xy[#inv_xy + 1] = p_rotate(th * i / z, p_rotate(angel, mirror(involute_t(r_TIF, ((AngelE) * j / steps))))) -- Mirroring
        end
    end
    inv_xy[#inv_xy + 1] = inv_xy[1] --table of values
    return linear_extrude(v(0, 0, b), inv_xy)
end

---- Function for the circle
function circle(r)
    local x, y = 0, 0
    local XY = {}
    for i = 1, 360 do
        local angle = i * math.pi / 180
        XY[i] = v(x + r * math.cos(angle), y + r * math.sin(angle))
    end
    return XY
end

---- For Internal Gear Extrude
function extrude(Contour, angle, dir_v, scale_v, z_steps)
    -- extrude a Contour to a shape by turning the contour to angle in z_steps
    -- extrude a Contour in a dircetion given by the vector dir_v
    -- extrude a Contour with a scaling factor given by vector scale_v
    -- Contour: a table of vectors as a closed contour (start point and end point is the same)
    -- angle: roation angle of contour along z_steps in deg
    -- dir_v: vector(x,y,z) direction of extrusion
    -- sacle_v: vector(x,y,z) scaling factors for extrudion
    -- z_steps: number of steps for the shape, mostly z_steps=2 if angle is equal zero

    local n_cont = #Contour
    local angle = angle / 180 * math.pi
    local Vertex = {}

    for j = 0, z_steps - 1 do
        local phi = angle * j / (z_steps - 1)
        local dir_vh = dir_v * j / (z_steps - 1)
        local scale_vh = (scale_v - v(1, 1, 1)) * (j / (z_steps - 1)) + v(1, 1, 1)
        for i = 1, n_cont - 1 do
            Vertex[i + j * n_cont] = v((dir_vh.x + scale_vh.x * (Contour[i].x * math.cos(phi) - Contour[i].y * math.sin(phi))),
                (dir_vh.y + scale_vh.y * (Contour[i].x * math.sin(phi) + Contour[i].y * math.cos(phi))),
                (dir_vh.z * scale_vh.z))
        end
        table.insert(Vertex, Vertex[1 + j * n_cont])
    end

    local vertex_sum_1 = v(0, 0, 0)
    local vertex_sum_m = v(0, 0, 0)

    for i = 1, n_cont - 1 do
        vertex_sum_1 = vertex_sum_1 + Vertex[i]
        vertex_sum_m = vertex_sum_m + Vertex[i + n_cont * (z_steps - 1)]
    end

    table.insert(Vertex, vertex_sum_1 / (n_cont - 1)) --n_cont*m_cont + 1
    table.insert(Vertex, vertex_sum_m / (n_cont - 1)) --n_cont*m_cont + 2

    Tri = {}
    local k = 1
    for j = 0, z_steps - 2 do
        for i = 0, n_cont - 2 do
            Tri[k] = v(i, i + 1, i + n_cont) + v(1, 1, 1) * n_cont * j
            Tri[k + 1] = v(i + 1, i + n_cont + 1, i + n_cont) + v(1, 1, 1) * n_cont * j
            k = k + 2
        end
    end
    for i = 0, n_cont - 2 do
        Tri[k] = v(i + 1, i, n_cont * z_steps)
        k = k + 1
    end
    for i = 0, n_cont - 2 do
        Tri[k] = v(i + n_cont * (z_steps - 1), i + 1 + n_cont * (z_steps - 1), n_cont * z_steps + 1)
        k = k + 1
    end
    return (polyhedron(Vertex, Tri))
end

---- Function Working Pressure angle
function WorkingAngelF(x) -- Acoording to 1992 [Harry Cheng], derivation of an explicit solution of the inverse involute function
    return (((math.pow(3 * x, (1 / 3)))) - (2 * x / 5) + (math.pow(9 / 175 * 3, (2 / 3))) * (math.pow(x, (5 / 3))) - (math.pow(2 / 175 * 3, (1 / 3))) * (math.pow(x, (7 / 3))) - ((144 / 67375) * (math.pow(x, (9 / 3))) + (3258 / 3128215) * (math.pow(3, (2 / 3))) * (math.pow(x, (11 / 3))) - (49711 / 153278125) * (math.pow(3, (1 / 3))) * (math.pow(x, (13 / 3))) - (1130112 / 9306171875) * (math.pow(x, (15 / 3))) + (5169659643 / 95304506171875) * (math.pow(3, (2 / 3))) * (math.pow(x, (17 / 3)))))
end

function center() --Calculation for the centre using the profile shift coefficients
    local alpha_rad = alpha_t * math.pi / 180;  -- Preassure angle
    local inv_a = math.tan(alpha_rad) - alpha_rad;  -- Involute function
    local inv_aw = ((2 * math.tan(alpha_rad) * (x_coef_ext - x_coef_int)) / (z_2 - z_1)) + inv_a;  -- Involute function working pressure angle

    -- Working preassure angle
    local alpha_aw = WorkingAngelF(inv_aw);
    -- Centre distance coeficiant factor
    local y = ((z_2 - z_1) * (math.cos(alpha_rad) - math.cos(alpha_aw))) / (2 * math.cos(alpha_aw));
    -- Center distance between Gears
    local a_x = (((z_2 - z_1) / 2) + y) * m;  -- Center distance
    meshDistance = a_x;
end

center();

function gear_formation()
    -- Formation of rotor gear
    local rotor_gear = gearProfile(z_2, m, alpha_t, x_coef_int, h_a_coef_p, h_f_coef_p, b);  -- Rotor gear
    local rot_rotor = rotate(0, 0, rotation);  -- Rotation of rotor gear
    emit(translate(0, 0, -0.1) * cylinder(r_a, 2), 1); -- Rotor gear Base formation
    emit(rot_rotor * difference(extrude(circle(r_a - 0.1), 0, v(0, 0, b), v(1, 1, 1), 20), rotor_gear), 1); -- Rotor gear formation
    set_brush_color(1, 0.3, 0.3, 0.4); -- Set color for rotor gear
	
    -- Parameters for the crescent calculation
    local crescent_outer_radius = r_TIF;  -- Radius of the rotor gear
    local crescent_root_radius = r_f;  -- Root radius of the rotor gear
	
    -- Formation of Gear housing
    local gear_housing = difference(cylinder(r_a + 16, b), cylinder(r_a + 1.75, b));  -- Gear housing cylinder
    local h = (z_2 * m) / 2 + m * 2 + 16;
    local inlet_trans = rotate(45, Z) * translate(h, 0, b / 2 + 0.2) * rotate(270, Y); -- Inlet translate
    local outlet_trans = rotate(45, Z) * translate(0, h, b / 2 + 0.2) * rotate(270, -X) * rotate(90, -Z); -- Outlet translate
    local incyl1 = inlet_trans * cylinder(b / 2 - 1, h);  -- Inlet cylinder 1
    local incyl2 = inlet_trans * translate(0, b - 2, 0) * cylinder(b / 2 - 1, h);  -- Inlet cylinder 2
    local incube = inlet_trans * translate(0, b / 2 - 1, 0) * cube(b - 2, b - 1, h); -- Inlet cube
    local inlet = union(union(incyl1,incyl2),incube); -- Inlet hole
    local outcyl1 = outlet_trans * cylinder(b / 2 - 1, h);  -- Outlet cylinder 1
    local outcyl2 = outlet_trans * translate(0, b - 2, 0) * cylinder(b / 2 - 1, h);  -- Outlet cylinder 2
    local outcube = outlet_trans * translate(0, b / 2 - 1, 0) * cube(b - 2, b - 1, h); -- Outlet cube
    local outlet = union(union(outcyl1,outcyl2),outcube); -- Outlet hole
    emit(intersection(difference(gear_housing, inlet), difference(gear_housing, outlet)), 0); -- Gear housing formation
    emit(translate(0, 0, -0.5) * cylinder(r_a + 16, 0.5), 0); -- Base for Gear housing formation
	
    -- Formation of idler gear
    local idler_gear = gearProfile(z_1, m, alpha_t, x_coef_ext, h_a_coef_p, h_f_coef_p, b);  -- Idler gear
    local rot_idler = rotate(0, 0, rotation * z_2 / z_1);  -- Rotation of idler gear
    emit(translate(0, meshDistance, 0) * rot_idler * cylinder(2 + r_b / 2, r_b + 15), 2); -- Shaft formation
    emit(translate(0, meshDistance, 0) * rot_idler * difference(idler_gear, extrude(circle(r_b / 2), 0, v(0, 0, b), v(1, 1, 1), 20)), 2); -- Idler gear formation
    set_brush_color(2, 0, 0, 0); -- Set color for idler gear
	
    -- Formation of the crescent
    local crescent_circle = translate(0, meshDistance, 0) * ccylinder(r_a + c * 2, b);  -- Cylinder for crescent
    local crescent_cube_right = translate(r_a + c, meshDistance + m * 2, 0) * cube(crescent_root_radius + c, crescent_root_radius + c, b + 0.1);  -- Right cube to remove sharp edges of the crescent
    local crescent_cube_left = translate(-(r_a), meshDistance + m * 2, 0) * cube(crescent_root_radius + c, crescent_root_radius + c, b + 0.1);  -- Left cube to remove sharp edges of the crescent
    local crescent_main = translate(0, 0, b / 2 + 0.1) * difference(ccylinder(crescent_outer_radius - c * 2, b), crescent_circle);  -- Crescent with sharp edges
    emit(intersection(difference(crescent_main, crescent_cube_right), difference(crescent_main, crescent_cube_left)), 0); -- Crescent formation without sharp edges
    set_brush_color(0, 0.2, 0.2, 0.2); -- Set color for crescent and gear housing
end

-- Gear Formation

gear_formation();
