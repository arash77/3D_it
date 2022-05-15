-------------------------------------------------------Internal Gear Pump--------------------------------------------------------

------------------------------------------------Declaration and Calculation End--------------------------------------------------
-- enable_variable_cache = true;                                                                                               -- caches the pointer

----------------------------------------------Passing parameters to user interface-----------------------------------------------
z_n=ui_numberBox("Number Of Teeth",20); 						                                                            -- Input for number of teeth, ideal: from 25
m=ui_numberBox("Module Of Gear",6);								                                                            -- Module of gear, Ideal Input: 3
alpha_t=ui_scalarBox("Pressure Angle(Deg)",20,1);				                                                            -- Pressure angle (changes the meshing parameters)
width=ui_numberBox("Width(mm)",10);								                                                            -- Width/Thickness of the gear
x_coef_int=ui_scalarBox("Internal Profile Shift(mm)",-0.2,0.1);	                                                            -- Profile shift factor for internal gear
x_coef_ext=ui_scalarBox("External Profile Shift(mm)",0,0.1);	                                                            -- Profile shift factor for external gear
h_a_coef_p=ui_scalarBox("Addendum Coefficient(mm)",1,0);					                                                -- Addendum height factor
h_f_coef_p=ui_scalarBox("Dedendum Coefficient(mm)",1.25,0);					                                                -- Dedendum height factor
rotation = ui_numberBox("Rotation",0);                                                                              		-- For rotation

--------------------------Function definition for Parametric equations for the involute profile points---------------------------
function tooth_involute(base_radius, inv_alpha)                                                                             -- inv_alpha = Involute angle
    return v(base_radius*(math.sin(inv_alpha) - inv_alpha*math.cos(inv_alpha)), 
			base_radius*(math.cos(inv_alpha) + inv_alpha* math.sin(inv_alpha)))
end

-----------------------------------Function defining for mirroring or inverting profile points-----------------------------------
function tooth_mirror(coord) 
	return v(-coord.x, coord.y) 
end

---------------------------------------------Function defining the rotational matrix---------------------------------------------
function rotate_points(angle, coord)                                                                                        -- angle = angle, coord = co-ordinates
    return v(math.cos(angle) * coord.x + math.sin(angle) * coord.y, math.cos(angle) * coord.y - math.sin(angle) * coord.x)
end

--------------------------------------Function defining angle between corresponding points---------------------------------------
function involute_angle(r_p1, r_p2)                                                                                       
	return (math.sqrt((r_p2 * r_p2 - r_p1 * r_p1) / (r_p1 * r_p1))) 
end

----------------------------Function Working Pressure angle -------------------
function wkp(x)
	return (((math.pow(3,(1/3)))*(math.pow(x,(1/3)))) - (2*x/5) + (math.pow(9/175*3,(2/3)))*(math.pow(x,(5/3))) - (math.pow(2/175*3, (1/3)))
						*(math.pow(x,(7/3))) - ((144/67375)*(math.pow(x,(9/3))) + (3258/3128215)*(math.pow(3,(2/3)))*(math.pow(x,(11/3)))
						- (49711/153278125)*(math.pow(3,(1/3)))*(math.pow(x,(13/3))) - (1130112/9306171875)*(math.pow(x,(15/3)))
						+ (5169659643/95304506171875)*(math.pow(3,(2/3)))*(math.pow(x,(17/3)))))
end

------------------------Function for the calculation for the required Internal and External Gear Profiles------------------------
function gear_profile(z, m, alpha_t, x_coef, h_a_coef, h_f_coef, width)
    local inv_xy = {}                                                                                                       -- Definition of the input parameters and calculation of the other base parameters for the involute profile.
    m_t = m;  						                                                                                        -- Module of gear
    alpha_t_rad = alpha_t*math.pi/180;                                                                                      -- Pressure angle
    z_t = z;  						                                                                                        -- Number of teeth
    x_coef = x_coef;				                                                                                        -- Profile shift co-efficient
    c = 0.167 * m_t;  			                                                                                            -- Clearance of tooth
	h_acoef = h_a_coef; 	                                                                                                -- Addendum coefficient
	h_dcoef = h_f_coef; 		                                                                                            -- Dedendum coefficient

-- Base Formula calculations required for the gear profile

-- Pitch Diameter:
    d_p = z_t * m_t 	   			                                                                                        -- Pitch Diamter
    r_p = d_p / 2 				   	                                                                                        -- Pitch radius

-- Base Diameter:
    d_b = d_p * math.cos(alpha_t_rad);                                                                                      -- Base diameter of gear
    r_b = d_b / 2 				                                                                                            -- Base radius

-- Addendum and Dedendum
	h_a = m_t* h_acoef; 			                                                                                        -- Addendum
	h_f = m_t* h_dcoef  			                                                                                        -- Dedendum

-- Tip Diameter / Addendum Diameter 
	d_a = m_t*(z_t+2*x_coef+2);                                                                                             -- Addendum diameter
    r_a = d_a / 2 				                                                                                            -- Addendum radius

-- Root Diameter / Dedendum Diameter
	d_f = m_t*z_t+ 2*x_coef*m_t - 2*(h_f)                                                                                   -- Root diamter
    r_f = d_f / 2 			                                                                                                -- Root_radius

-- True Involute Diameter:
    Q_fp = ((m_t*((math.pi/4)-(math.tan(alpha_t_rad)))-c*math.tan(alpha_t_rad))*(1+math.sin(alpha_t_rad)))/(math.cos(alpha_t_rad))
	d_TIF = math.sqrt(math.pow(d_p * math.sin(alpha_t_rad) - 2 *(h_a - (m_t * x_coef) - Q_fp *(1 - math.sin(alpha_t_rad)))/(math.sin(alpha_t_rad)), 2) +   d_b * d_b)    -- True involute diameter
    r_TIF = d_TIF / 2; 				                                                                                        -- true form radius 

-- Tooth thickness on the pitch Circle
    S_0 = m_t*((math.pi/2)+ 2*x_coef* math.tan(alpha_t_rad));

-- Involute function
    inv_a = math.tan(alpha_t_rad) - alpha_t_rad; 

-- Tooth thickness on the form circle
    alpha_f = math.acos((d_p*math.cos(alpha_t_rad))/d_TIF); 				                                                -- Involute angle along form circle
    inv_Ff  = math.tan(alpha_f) - alpha_f;									                                                -- Involute function
    s_f     = d_TIF*((S_0/d_p)+inv_a - inv_Ff)+(x_coef*math.tan(alpha_f));                                                  -- Thickness of the gear teeth along form circle 
    omega   = (s_f)/(0.5*d_TIF); 											                                                -- Angle swept along form circle for corresponding thickeness

-- Function for stariing and ending of involute between two radius
    tooth_ang = (((math.pi * m_t / 2) + 2 * m_t * x_coef * math.tan(alpha_t_rad)) /
                    r_p + 2 * math.tan(alpha_t_rad) - 2 * alpha_t_rad)
    res = 30;                                                                                                               -- Defining iterating points 

-- Iterating the points for the involute points and iterating teeth around the circle
    for i = 1, z_t do
        th = 2 * math.pi 				                                                                                    -- Iteration angle of the teeth around the gear diameter
        th1 = involute_angle(r_b, r_TIF);                                                                                   -- Start angle to define the start of the involute curve 
        th2 = involute_angle(r_TIF, r_a);                                                                                   -- End asngle to determine the end of the involute curve 

-- Defining iterations for the involute teeth profile 
        for j = 1, res do
			inv_xy[#inv_xy + 1] = rotate_points(th*i/z_t, tooth_involute(r_TIF, ((th2) * j / res)))                         -- The one side of the involute points are obtained 
        end
        for j = res, 1, -1 do
            inv_xy[#inv_xy + 1] = rotate_points(th*i/z_t, rotate_points(omega, tooth_mirror(tooth_involute(r_TIF, ((th2)* j / res)))))       -- The second side of the involute points are obtained using mirroring function 
        end
    end
    inv_xy[#inv_xy + 1] = inv_xy[1] 	                                                                                    -- Used to generate a table of values 
	interior = linear_extrude(v(0,0,width),inv_xy)
    return interior
end
------------------------------Defining the Default Parameters for the Internal Gear Pump Components------------------------------
function gear(params)
      return gear_profile(
			params.z or 20,				                                                                                    -- Number of teeth of the Internal Gear
			params.m or 3, 			                                                                                        -- Module  
			params.alpha_t or 28, 		                                                                                    -- Pressure Angle
            params.x_coef or 0, 		                                                                                    -- Profile shift coefficient 
			params.h_a_coef or 1,		                                                                                    -- Addendum coefficient
			params.h_f_coef or 1.25,                                                                                        -- Dedendum coefficient
			params.width or 15)			                                                                                    -- Thickness/Width 
end

-----------------------------------Function defining the parametrical equation for the circle------------------------------------
function circle(r)                                                                                                          -- r -radius of circle
  local x, y = 0, 0
  local XY={}
  for i = 1, 360 do
    local angle = i * math.pi / 180
    XY[i] = v(x + r * math.cos( angle ), y + r * math.sin( angle ))
  end
  return XY
end

----------------------------------------------------For Internal Gear Extrude----------------------------------------------------
function extrude(Contour, angle, dir_v, scale_v, z_steps)
-- extrude a Contour to a shape by turning the contour to angle in z_steps
-- extrude a Contour in a dircetion given by the vector dir_v
-- extrude a Contour with a scaling factor given by vector scale_v 
-- Contour: a table of vectors as a closed contour (start point and end point is the same)
-- angle: roation angle of contour along z_steps in deg
-- dir_v: vector(x,y,z) direction of extrusion
-- sacle_v: vector(x,y,z) scaling factors for extrudion
-- z_steps: number of steps for the shape, mostly z_steps=2 if angle is equal zero

   local n_cont= #Contour
   local angle= angle/180*math.pi
   local Vertex= {}

   for j= 0,z_steps-1 do
      local phi= angle*j/(z_steps-1)
      local dir_vh= dir_v*j/(z_steps-1)
      local scale_vh= (scale_v - v(1,1,1))*(j/(z_steps-1)) + v(1,1,1)
	  for i= 1,n_cont-1 do
          Vertex[i+j*n_cont]= v((dir_vh.x + scale_vh.x * (Contour[i].x*math.cos(phi) - Contour[i].y*math.sin(phi))),
                                (dir_vh.y + scale_vh.y * (Contour[i].x*math.sin(phi) + Contour[i].y*math.cos(phi))),
                                (dir_vh.z * scale_vh.z))
      end
      table.insert(Vertex,Vertex[1+j*n_cont])
   end

   local vertex_sum_1 = v(0,0,0)
   local vertex_sum_m = v(0,0,0)

   for i= 1,n_cont-1 do
      vertex_sum_1= vertex_sum_1 + Vertex[i]
      vertex_sum_m= vertex_sum_m + Vertex[i+n_cont*(z_steps-1)]
   end

   table.insert(Vertex,vertex_sum_1/(n_cont-1))                                                                             --n_cont*m_cont + 1
   table.insert(Vertex,vertex_sum_m/(n_cont-1))                                                                             --n_cont*m_cont + 2

   Tri= {}                                                                                                                  -- !!! the index on table with Vertex starts with zero !!!
   local k= 1

   for j=0,z_steps-2 do
      for i= 0,n_cont-2 do
         Tri[k]=   v(i, i+1, i+n_cont) + v(1,1,1)*n_cont*j
         Tri[k+1]= v(i+1, i+n_cont+1, i+n_cont) + v(1,1,1)*n_cont*j
         k= k+2
      end
   end
   for i= 0,n_cont-2 do
      Tri[k]= v(i+1,i,n_cont*z_steps)
      k= k+1
   end
   for i= 0,n_cont-2 do
      Tri[k]= v(i+n_cont*(z_steps-1),i+1+n_cont*(z_steps-1),n_cont*z_steps+1)
      k= k+1
   end
   return(polyhedron(Vertex,Tri))
end

x_coef_ext = x_coef_int + 0.1;
if alpha_t < 4 then
	alpha_t = 5;
end

----------------------------------------------------External Gear Formation------------------------------------------------------
externalGear = gear({z=z_n;m=m;alpha_t=alpha_t;x_coef=x_coef_int;h_a_coef=h_a_coef_p;h_f_coef=h_f_coef_p;width=width})

---------------------Calculation for the centre distance between gears using the profile shift coefficients----------------------
z_2 = z_n;                 		                                                                                            -- Number of teeth Internal Gear
z_1 = z_n-8;				 		                                                                                        -- Number of teeth External Gear
alpha_rad = alpha_t*math.pi/180                                                                                             -- Preassure angle
inv_a = math.tan(alpha_rad) - alpha_rad;                                                                                    -- Involute function
inv_aw = ((2*math.tan(alpha_rad) * (x_coef_int - x_coef_ext))/(z_2 - z_1)) + inv_a;                                         -- Involute function working pressure angle
alpha_aw = wkp(inv_aw)                                                                                                      -- Working pressure angle
y_c = ((z_2 - z_1) * (math.cos(alpha_rad) - math.cos(alpha_aw)))/ (2* math.cos(alpha_aw));                                  -- Centre distance incremental factor
meshDistance = (((z_2 - z_1)/2)+y_c) * m;				                                                                    -- Center distance is assagined for cresant calculation
addOn_distance= m * 2;
crescent_outerRadius=r_TIF;		                                                                                            -- Radius of the internal gear for the cylinder 
crescent_clearanceOut=c;		                                                                                            -- Clearance of the internal gear for the cylinder
crescent_rootRadius=r_f;		                                                                                            -- Root radius of the internal gear
internal_radius=r_a;																										-- Pitch radius of internal gear
crescent_cylTop=ccylinder(crescent_outerRadius-crescent_clearanceOut*2,width);	                                            -- Top cylinder for crescent 
radius=(z_n*m)/2;

----------------------------------------------------Internal Gear Formation------------------------------------------------------
internalGear = gear({z=z_n-8;m=m;alpha_t=alpha_t;x_coef=x_coef_ext;h_a_coef=h_a_coef_p;h_f_coef=h_f_coef_p;width=width})

--------------------------------------------------------Crescent Filler----------------------------------------------------------
crescent_innerRadius=r_a;   			                                                                                    -- Radius of the external gear for the cylinder
crescent_clearanceIn=c;					                                                                                    -- Clearance of the external gear for the cylinder
crescent_translation=meshDistance;		                                                                                    -- Translation for crescent formation
crescent_cylBottom=translate(0,crescent_translation,0)*ccylinder(crescent_innerRadius+crescent_clearanceIn*2,width); 	    --Inner circle for internal gear formation
crescent_cubeRight=translate(crescent_innerRadius+crescent_clearanceOut,meshDistance+m*2,0)*
                   cube(crescent_rootRadius+crescent_clearanceOut, crescent_rootRadius+crescent_clearanceOut,width+0.1)     -- Cube to remove sharp edges of the crescent
crescent_cubeLeft=translate(-(crescent_innerRadius+crescent_clearanceOut),meshDistance+m*2,0)*
                   cube(crescent_rootRadius+crescent_clearanceOut, crescent_rootRadius+crescent_clearanceOut,width+0.1)     -- Cube to remove sharp edges of the crescent
crescent_Main=translate(0,0,width/2+0.1)*difference(crescent_cylTop,crescent_cylBottom)                                     -- Crescent Filler Formation with sharp edges
crescent_Full=intersection(difference(crescent_Main,crescent_cubeRight),difference(crescent_Main,crescent_cubeLeft))        -- Crescent Filler Formation without sharp edges
emit(crescent_Full,0)                                                                                                         -- Crescent Filler Formation
set_brush_color(0,0.2,0.2,0.2)

-----------------------------------Emitting the Internal Gear, the External Gear, and shaft--------------------------------------



outer_circle= extrude(circle(internal_radius-0.1), 0, v(0,0,width), v(1,1,1), 20)                                     		-- Outer circle for External gear formation
r1 = rotate(0,0, rotation)                                                                                                  -- External Gear Rotation
emit(r1*difference(outer_circle,externalGear),1);                                                                           -- External gear formation
emit(translate(0,0,-0.1)*cylinder(internal_radius,2),1)															    		-- External gear Base formation
set_brush_color(1,0.1,0.2,0.2)

flangeCube=translate(0,-meshDistance,0)*cube(addOn_distance*2,addOn_distance*2,width)
flangeCylinder=cylinder(addOn_distance*2,width+20)
cube_flange=translate(0,-meshDistance,0)*cube(addOn_distance*2,addOn_distance*2,width+5)
flangeDiff=difference(extrude(circle(addOn_distance*2), 0, v(0,0,width), v(1,1,1), 20),flangeCube)
flangeFull=difference(internalGear,flangeDiff)
r2 = rotate(0,0, rotation*z_2/z_1)                                                                                          -- Internal Gear Rotation
emit(translate(0,meshDistance,0.1)*r2*difference(flangeCylinder,cube_flange),2)									    	 	-- Shaft Formation
emit(translate(0,meshDistance,0)*r2*difference(flangeFull,flangeDiff),2)                         						 	-- Internal Gear Formation
set_brush_color(2,0,0,0)


