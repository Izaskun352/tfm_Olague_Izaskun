
# pyMOL

#dnah11
#input proteina - le llamamos 'prot'
load /Users/Usuario/Desktop/IZASKUN/TFM/Output/Gráficos/Alphafold/Q96DT5_38_4516_8j07.908.A.cif , prot; bg_color white

color wheat, prot;

#seleccionamos dominios
select dominio_stem, prot and resi 1-1854;
color forest, dominio_stem;

select dominio_AAA1, prot and resi 1274-2076;
color tv_blue, dominio_AAA1;

select dominio_AAA2, prot and resi 2136-2366;
color marine, dominio_AAA2;

select dominio_AAA3, prot and resi 2472-2719;
color slate, dominio_AAA3;

select dominio_AAA4, prot and resi 2817-3066;
color lightblue, dominio_AAA4;

select dominio_stalk, prot and resi 3072-3403;
color chartreuse, dominio_stalk;

select dominio_AAA5, prot and resi 3459-3686;
color skyblue, dominio_AAA5;

select dominio_AAA6, prot and resi 3896-4122;
color deepblue, dominio_AAA6;

select coiled_coil1, prot and resi 1274-1327;
color violet, coiled_coil1;

select coiled_coil2, prot and resi 3072-3136;
color violet, coiled_coil2;

select coiled_coil3, prot and resi 3312-3403;
color violet, coiled_coil3;

select coiled_coil4, prot and resi 3668-3703;
color violet, coiled_coil4;

#seleccionamos variantes
select variant_c1, prot and resi 1998+2002+2006;
color tv_red, variant_c1;

select variant_c2, prot and resi 2060+2065+2077;
color tv_red, variant_c2;

select variant_c3, prot and resi 1939+1946+1949;
color tv_red, variant_c3;

select variant_c4, prot and resi 4141+4154;
color tv_red, variant_c4;

select variant_c5, prot and resi 2638+2639;
color tv_red, variant_c5;

select variant_stem, prot and resi 723+891+1343+1610;
color tv_red, variant_stem;

select variant_stalk, prot and resi 3162+3385;
color tv_red, variant_stalk;

select variant_c6, prot and resi 3523+3650+3385+4435;
color tv_red, variant_c6;

show surface, prot
show licorice, variant_c1
show spheres, variant_c1;
show spheres, variant_c2;
show spheres, variant_c3;
show spheres, variant_c4;
show spheres, variant_stem;
show spheres, variant_c5;
show spheres, variant_c6

show sticks, variant_c1; set stick_ball, on; set stick_ball_ratio, 1.0

set transparency, 0.7
set transparency, 0.6, polymer.protein

#exportar
set ray_opaque_background, off; png /Users/Usuario/Desktop/IZASKUN/TFM/Output/pymol/dnah11_full.png, width=2500, height=2000, dpi=300, ray=1

#complejo en pulmo de humanos
#input proteina - le llamamos 'complex'
load /Users/Usuario/Desktop/IZASKUN/TFM/Output/Gráficos/Alphafold/8j07_updated.cif , complex; bg_color white
color darksalmon, complex;

#coloreamos cadena
select dnah9_1, chain n9
color tv_blue, dnah9_1;

select dnah9_2, chain p9
color marine, dnah9_2;

select dnah9_3, chain r9;
color slate, dnah9_3

select dnah9_4, chain t9;
color lightblue, dnah9_4


#coloreamos variantes en la cadena
color red, dnah9_1 and resi 3761;
show sticks, dnah9_1 and resi 3761

color red, dnah9_1 and resi 3350;
show sticks, dnah9_1 and resi 3350

color red, dnah9_1 and resi 1842;
show sticks, dnah9_1 and resi 1842

color red, dnah9_1 and resi 2072;
show sticks, dnah9_1 and resi 2072

color red, dnah9_1 and resi 2156;
show sticks, dnah9_1 and resi 2156

color red, dnah9_1 and resi 1198;
show sticks, dnah9_1 and resi 1198

color orange, dnah9_1 and resi 4141;
show sticks, dnah9_1 and resi 4141

