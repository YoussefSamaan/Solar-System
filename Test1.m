clc
clear

%Planet's masses in Kg
Sun_mass = 1.989*10^32;
Mercury_mass = 3.3*10^23;
Venus_mass = 4.87*10^24; 
Earth_mass = 5.97*10^24;
Mars_mass = 6.42*10^23;
Jupiter_mass = 1.90*10^27; 
Saturn_mass = 5.68*10^26; 
Uranus_mass = 8.68*10^25; 
Neptune_mass = 1.02*0^26;

%Planet's average distance (radius from the sun) in meters (initially on the x-axis (y=0))
Sun_radius = 0;
Mercury_radius = 5.79*10^10;
Venus_radius = 1.08*10^11;
Earth_radius = 1.496*10^11;
Mars_radius = 2.28*10^11;
Jupiter_radius = 7.78*10^11;
Saturn_radius = 1.43*10^12;
Uranus_radius = 2.87*10^12;
Neptune_radius = 4.50*10^12;

%Planet's velocity i meters/seconds
Sun_Vi = 0;
Mercury_Vi = 47900;
Venus_Vi = 35000;
Earth_Vi = 29800;
Mars_Vi = 24100;
Jupiter_Vi = 30100;
Saturn_Vi = 40700;
Uranus_Vi = 50800;
Neptune_Vi = 90400;

Planets_radii = [6.96*10^8,
                2.44*10^6,
                6.05*10^6,
                6.37*10^6,
                3.39*10^6,
                6.99*10^7,
                5.82*10^7,
                2.54*10^7,
                2.46*10^7];

save 'Test4.mat'