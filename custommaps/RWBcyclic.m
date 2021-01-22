function cm_data=RWBcyclic(m)
cm = [0	0	0
0	0.0298039223998785	0.0494117662310600
0	0.0596078447997570	0.0988235324621201
0	0.0894117653369904	0.148235291242600
0	0.119215689599514	0.197647064924240
0	0.149019613862038	0.247058823704720
0	0.178823530673981	0.296470582485199
0	0.208627447485924	0.345882356166840
0	0.238431379199028	0.395294129848480
0	0.268235296010971	0.444705903530121
0	0.298039227724075	0.494117647409439
0	0.327843129634857	0.543529450893402
0	0.357647061347961	0.592941164970398
0	0.387450993061066	0.642352938652039
0	0.417254894971848	0.691764712333679
0	0.447058826684952	0.741176486015320
0.0625000000000000	0.481617659330368	0.757352948188782
0.125000000000000	0.516176462173462	0.773529410362244
0.187500000000000	0.550735294818878	0.789705872535706
0.250000000000000	0.585294127464294	0.805882334709168
0.312500000000000	0.619852960109711	0.822058856487274
0.375000000000000	0.654411792755127	0.838235318660736
0.437500000000000	0.688970565795898	0.854411780834198
0.500000000000000	0.723529398441315	0.870588243007660
0.562500000000000	0.758088231086731	0.886764705181122
0.625000000000000	0.792647063732147	0.902941167354584
0.687500000000000	0.827205896377564	0.919117629528046
0.750000000000000	0.861764729022980	0.935294151306152
0.812500000000000	0.896323502063751	0.951470613479614
0.875000000000000	0.930882334709168	0.967647075653076
0.937500000000000	0.965441167354584	0.983823537826538
1	1	1
0.990686297416687	0.957843124866486	0.943627476692200
0.981372535228729	0.915686249732971	0.887254893779755
0.972058832645416	0.873529434204102	0.830882370471954
0.962745070457459	0.831372559070587	0.774509787559509
0.953431367874146	0.789215683937073	0.718137264251709
0.944117665290833	0.747058808803558	0.661764681339264
0.934803903102875	0.704901993274689	0.605392158031464
0.925490200519562	0.662745118141174	0.549019634723663
0.916176497936249	0.620588243007660	0.492647051811218
0.906862735748291	0.578431367874146	0.436274498701096
0.897549033164978	0.536274492740631	0.379901975393295
0.888235330581665	0.494117647409439	0.323529422283173
0.878921568393707	0.451960802078247	0.267156869173050
0.869607865810394	0.409803926944733	0.210784316062927
0.860294103622437	0.367647081613541	0.154411762952805
0.850980401039124	0.325490206480026	0.0980392172932625
0.797794103622437	0.305147081613541	0.0919117629528046
0.744607865810394	0.284803926944733	0.0857843160629273
0.691421568393707	0.264460802078247	0.0796568617224693
0.638235330581665	0.244117647409439	0.0735294148325920
0.585049033164978	0.223774522542954	0.0674019604921341
0.531862735748291	0.203431382775307	0.0612745098769665
0.478676468133926	0.183088243007660	0.0551470592617989
0.425490200519562	0.162745103240013	0.0490196086466312
0.372303932905197	0.142401963472366	0.0428921580314636
0.319117665290833	0.122058823704720	0.0367647074162960
0.265931367874146	0.101715691387653	0.0306372549384832
0.212745100259781	0.0813725516200066	0.0245098043233156
0.159558832645416	0.0610294118523598	0.0183823537081480
0.106372550129890	0.0406862758100033	0.0122549021616578
0.0531862750649452	0.0203431379050016	0.00612745108082891
0	0	0];

if nargin < 1
    cm_data = cm;
else
    hsv=rgb2hsv(cm);
    cm_data=interp1(linspace(0,1,size(cm,1)),hsv,linspace(0,1,m));
    cm_data=hsv2rgb(cm_data);
  
end
end