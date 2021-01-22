function cm_data=coldintensity(m)
cm = [[0,0,0;0.000392156885936856,0.0325490199029446,0.0439215674996376;0.000784313771873713,0.0650980398058891,0.0878431349992752;0.00117647065781057,0.0976470634341240,0.131764709949493;0.00156862754374743,0.130196079611778,0.175686269998550;0.00196078442968428,0.162745103240013,0.219607844948769;0.00235294131562114,0.195294126868248,0.263529419898987;0.00274509820155799,0.227843150496483,0.307450979948044;0.00313725508749485,0.260392159223557,0.351372539997101;0.00352941197343171,0.292941182851791,0.395294129848480;0.00392156885936856,0.325490206480026,0.439215689897537;0.00371517054736614,0.336842119693756,0.455108374357224;0.00350877223536372,0.348194032907486,0.471001029014587;0.00330237369053066,0.359545946121216,0.486893713474274;0.00309597537852824,0.370897859334946,0.502786397933960;0.00288957706652582,0.382249742746353,0.518679082393646;0.00268317875452340,0.393601655960083,0.534571707248688;0.00247678044252098,0.404953569173813,0.550464391708374;0.00227038189768791,0.416305482387543,0.566357076168060;0.00206398358568549,0.427657395601273,0.582249760627747;0.00185758527368307,0.439009308815002,0.598142445087433;0.00165118684526533,0.450361222028732,0.614035069942474;0.00144478853326291,0.461713135242462,0.629927754402161;0.00123839022126049,0.473065048456192,0.645820438861847;0.00103199179284275,0.484416961669922,0.661713123321533;0.000825593422632664,0.495768845081329,0.677605807781220;0.000619195110630244,0.507120788097382,0.693498492240906;0.000412796711316332,0.518472671508789,0.709391117095947;0.000206398355658166,0.529824614524841,0.725283801555634;0,0.541176497936249,0.741176486015320;0.0325259529054165,0.553633272647858,0.749942362308502;0.0650519058108330,0.566089987754822,0.758708178997040;0.0975778624415398,0.578546762466431,0.767474055290222;0.130103811621666,0.591003477573395,0.776239931583405;0.162629768252373,0.603460252285004,0.785005807876587;0.195155724883080,0.615916967391968,0.793771624565125;0.227681666612625,0.628373742103577,0.802537500858307;0.260207623243332,0.640830457210541,0.811303377151489;0.292733579874039,0.653287231922150,0.820069193840027;0.325259536504746,0.665743947029114,0.828835070133209;0.357785493135452,0.678200721740723,0.837600946426392;0.390311449766159,0.690657436847687,0.846366763114929;0.422837376594543,0.703114211559296,0.855132639408112;0.455363333225250,0.715570926666260,0.863898515701294;0.487889289855957,0.728027701377869,0.872664391994476;0.520415246486664,0.740484416484833,0.881430208683014;0.552941203117371,0.752941191196442,0.890196084976196;0.574625194072723,0.764705896377564,0.894809722900391;0.596309125423431,0.776470601558685,0.899423301219940;0.617993116378784,0.788235306739807,0.904036939144135;0.639677047729492,0.800000011920929,0.908650517463684;0.661361038684845,0.811764717102051,0.913264155387878;0.683045029640198,0.823529422283173,0.917877733707428;0.704728960990906,0.835294127464294,0.922491371631622;0.726412951946259,0.847058832645416,0.927104949951172;0.748096883296967,0.858823537826538,0.931718587875366;0.769780874252319,0.870588243007660,0.936332166194916;0.791464805603027,0.882352948188782,0.940945804119110;0.813148796558380,0.894117653369904,0.945559382438660;0.834832787513733,0.905882358551025,0.950173020362854;0.856516718864441,0.917647063732147,0.954786598682404;0.878200709819794,0.929411768913269,0.959400236606598;0.899884641170502,0.941176474094391,0.964013814926148;0.921568632125855,0.952941179275513,0.968627452850342]];

if nargin < 1
    cm_data = cm;
else
    hsv=rgb2hsv(cm);
    cm_data=interp1(linspace(0,1,size(cm,1)),hsv,linspace(0,1,m));
    cm_data=hsv2rgb(cm_data);
  
end
end