path_to_include= 'C:\Users\njgreen3\Documents\MATLAB\CoolPropStuff'; %specify path to coolprop's include folder

% Loading shared library
if ~libisloaded('coolprop') %checking whether library is already loaded
%     addpath(path_to_lib)
    addpath(path_to_include)
    libname = 'libCoolProp'; % OSX and linux
    if ispc
        libname = 'CoolProp';
    end
    loadlibrary(libname,'CoolPropLib.h','includepath',path_to_include,'alias','coolprop'); % loading library with alias coolprop
    disp('loaded CoolProp shared library.')
    disp('loaded these functions: ')
    libfunctions coolprop
end

% set error handling variables
buffer_size = 255;
ierr = 0;
b = (1:1:buffer_size);
herr= char(b);

% Close existing handles (to keep things tidy)
if exist('handle_source','var')
    calllib('coolprop','AbstractState_free', handle_source, ierr,herr,buffer_size);
end
if exist('handle_sink','var')
    calllib('coolprop','AbstractState_free', handle_sink, ierr,herr,buffer_size);
end
if exist('handle_wf','var')
    calllib('coolprop','AbstractState_free', handle_wf, ierr,herr,buffer_size);
end

% define the source, sink, and working fluids and set up respective handles
fluid_source = 'water';
fluid_sink = 'water';
fluid_wf = 'R245fa';
backend = 'HEOS';
[handle_source, ~] = calllib('coolprop','AbstractState_factory',backend,fluid_source,ierr,herr,buffer_size);
[handle_sink, ~] = calllib('coolprop','AbstractState_factory',backend,fluid_sink,ierr,herr,buffer_size);
[handle_wf, ~] = calllib('coolprop','AbstractState_factory',backend,fluid_wf,ierr,herr,buffer_size);

Tin_source = 81.3 + 273.15;
Tin_sink = 8.5 + 273.15;
m_dot_source = 937.44*.971/60;  %   937 l/min * .9710 kg/l (@ 81.3 C) * 1 min/ 60 s = 15 kg/s
m_dot_sink = 937.44*.9998/60;  %   937 l/min * .9998 kg/l (@ 8.5 C) * 1 min/ 60 s = 15 kg/s
m_dot_wf_init = 3;
m_dot_wf_max = 12;
H_wf_init = 2.3e+05;
p_hi = 1e6;
p_low = 1.5e5;
p_atm = 101325;

U_evap = 1500;
A_evap = 20.3;
U_cond = 1400;
A_cond = 20.6;

eff_turbine = 0.78;
eff_pump = 0.78;
eff_inverter = 0.93;
pf_grid = 0.9;

P_spd_speed = [834.591559026165	844.140855611678	853.690152197192	863.239448782706	872.788745368220	882.338041953733	891.887338539247	901.436635124761	910.985931710274	920.535228295788	930.084524881302	939.633821466816	949.183118052329	958.732414637843	968.281711223357	977.831007808871];
P_spd_power = -1/3*[166663.989372846	163937.628261044	155852.617467788	140628.948970298	116654.771019611	83123.6124212207	40806.0621676137	-7481.19463142690	-57239.4174183274	-103560.341265182	-142633.126036591	-172599.754334336	-193461.262910890	-206409.980979167	-213110.460435838	-215209.127126632];

magCurveX = [4.38180036558486,4.38174885826274,4.38160514998934,4.38136897475881,4.38104009678620,4.38061831012584,4.38010343828725,4.37949533384872,4.37879387806883,4.37799898049604,4.37711057857683,4.37612863726219,4.37505314861309,4.37388413140472,4.37262163073000,4.37126571760228,4.36981648855757,4.36827406525629,4.36663859408486,4.36491024575712,4.36308921491577,4.36117571973394,4.35917000151707,4.35707232430512,4.35488297447523,4.35260226034509,4.35023051177681,4.34776807978174,4.34521533612607,4.34257267293736,4.33984050231218,4.33701925592483,4.33410938463722,4.33111135811007,4.32802566441539,4.32485280965036,4.32159331755274,4.31824772911768,4.31481660221621,4.31130051121529,4.30770004659964,4.30401581459516,4.30024843679427,4.29639854978303,4.29246680477012,4.28845386721775,4.28436041647457,4.28018714541048,4.27593476005362,4.27160397922927,4.26719553420104,4.26271016831400,4.25814863664019,4.25351170562617,4.24880015274289,4.24401476613782,4.23915634428933,4.23422569566346,4.22922363837290,4.22415099983852,4.21900861645316,4.21379733324782,4.20851800356043,4.20317148870689,4.19775865765474,4.19228038669922,4.18673755914192,4.18113106497192,4.17546180054948,4.16973066829232,4.16393857636446,4.15808643836768,4.15217517303557,4.14620570393024,4.14017895914162,4.13409587098948,4.12795737572805,4.12176441325340,4.11551792681341,4.10921886272058,4.10286817006739,4.09646680044454,4.09001570766185,4.08351584747190,4.07696817729644,4.07037365595562,4.06373324339989,4.05704790044477,4.05031858850841,4.04354626935185,4.03673190482225,4.02987645659881,4.02298088594153,4.01604615344290,4.00907321878227,4.00206304048323,3.99501657567369,3.98793477984894,3.98081860663751,3.97366900756991,3.96648693185024,3.95927332613071,3.95202913428903,3.94475529720864,3.93745275256197,3.93012243459643,3.92276527392346,3.91538219731040,3.90797412747527,3.90054198288457,3.89308667755385,3.88560912085135,3.87811021730454,3.87059086640945,3.86305196244319,3.85549439427921,3.84791904520559,3.84032679274629,3.83271850848528,3.82509505789374,3.81745730016010,3.80980608802314,3.80214226760798,3.79446667826510,3.78678015241228,3.77908351537953,3.77137758525703,3.76366317274601,3.75594108101258,3.74821210554465,3.74047703401169,3.73273664612764,3.72499171351666,3.71724299958197,3.70949125937763,3.70173723948339,3.69398167788239,3.68622530384205,3.67846883779782,3.67071299123997,3.66295846660343,3.65520595716059,3.64745614691708,3.63970971051069,3.63196731311310,3.62422961033481,3.61649724813300,3.60877086272239,3.60105108048917,3.59333851790788,3.58563378146139,3.57793746756386,3.57025016248670,3.56257244228758,3.55490487274249,3.54724800928078,3.53960239692323,3.53196857022320,3.52434705321072,3.51673835933975,3.50914299143828,3.50156144166164,3.49399419144875,3.48644171148140,3.47890446164664,3.47138289100211,3.46387743774446,3.45638852918082,3.44891658170326,3.44146200076637,3.43402518086771,3.42660650553158,3.41920634729555,3.41182506770016,3.40446301728174,3.39712053556806,3.38979795107731,3.38249558131982,3.37521373280302,3.36795270103945,3.36071277055771,3.35349421491655,3.34629729672193,3.33912226764720,3.33196936845631,3.32483882903002,3.31773086839521,3.31064569475728,3.30358350553548,3.29654448740137,3.28952881632041,3.28253665759643,3.27556816591922,3.26862348541531,3.26170274970155,3.25480608194199,3.24793359490767,3.24108539103946,3.23426156251409,3.22746219131302,3.22068734929464,3.21393709826922,3.20721149007718,3.20051056667031,3.19383436019597,3.18718289308451,3.18055617813962,3.17395421863181,3.16737700839487,3.16082453192555,3.15429676448609,3.14779367220992,3.14131521221046,3.13486133269296,3.12843197306926,3.12202706407586,3.11564652789477,3.10929027827776,3.10295822067325,3.09665025235667,3.09036626256365,3.08410613262631,3.07786973611266,3.07165693896899,3.06546759966537,3.05930156934432,3.05315869197229,3.04703880449436,3.04094173699213,3.03486731284433,3.02881534889093,3.02278565559987,3.01677803723723,3.01079229204015,3.00482821239321,2.99888558500733,2.99296419110234,2.98706380659202,2.98118420227280,2.97532514401500,2.96948639295751,2.96366770570527,2.95786883453007,2.95208952757415,2.94632952905715,2.94058857948581,2.93486641586706,2.92916277192392,2.92347737831462,2.91780996285473,2.91216025074234,2.90652796478640,2.90091282563796,2.89531455202468,2.88973286098822,2.88416746812489,2.87861808782913,2.87308443354032,2.86756621799248,2.86206315346701,2.85657495204877,2.85110132588486,2.84564198744669,2.84019664979513,2.83476502684866,2.82934683365453,2.82394178666318,2.81854960400551,2.81317000577336,2.80780271430307,2.80244745446201,2.79710395393820,2.79177194353315,2.78645115745747,2.78114133362996,2.77584221397934,2.77055354474936,2.76527507680686,2.76000656595287,2.75474777323686,2.74949846527404,2.74425841456572,2.73902739982260,2.73380520629151,2.72859162608476,2.72338645851286,2.71818951042027,2.71300059652414,2.70781953975611,2.70264617160733,2.69748033247639,2.69232187202042,2.68717064950916,2.68202653418229,2.67688940560955,2.67175915405425,2.66663568083954,2.66151889871806,2.65640873224441,2.65130511815072,2.64620800572550,2.64111735719539,2.63603314810990,2.63095536772951,2.62588401941656,2.62081912102920,2.61576070531882,2.61070882033003,2.60566352980398,2.60062491358485,2.59559306802911,2.59056810641810,2.58555015937357,2.58053937527642,2.57553592068831,2.57053998077642,2.56555175974151,2.56057148124860,2.55559938886119,2.55063574647817,2.54568083877413,2.54073497164246,2.53579847264172,2.53087169144494,2.52595500029222,2.52104879444598,2.51615349264983,2.51126953759002,2.50639739636032,2.50153756092972,2.49669054861341,2.49185690254661,2.48703719216173,2.48223201366831,2.47744199053644,2.47266777398266,2.46791004345961,2.46316950714822,2.45844690245322,2.45374299650177,2.44905858664485,2.44439450096223,2.43975159876997,2.43513077113142,2.43053294137099,2.42595906559133,2.42141013319314,2.41688716739841,2.41239122577668,2.40792340077419,2.40348482024621,2.39907664799263,2.39470008429621,2.39035636646443,2.38604676937395,2.38177260601811,2.37753522805820,2.37333602637694,2.36917643163557,2.36505791483381,2.36098198787297,2.35695020412195,2.35296415898681,2.34902549048260,2.34513587980922,2.34129705192953,2.33751077615080,2.33377886670975,2.33010318335979,2.32648563196198,2.32292816507871,2.31943278257071,2.31600153219686,2.31263651021747,2.30933986200010,2.30611378262909,2.30296051751736,2.29988236302219,2.29688166706332,2.29396082974464,2.29112230397856,2.28836859611377,2.28570226656596,2.28312593045133,2.28064225822380];
magCurveV = [98.1532063054068,98.8036071709152,99.4518999785459,100.098037520928,100.741973127402,101.383660677656,102.023054614968,102.660109959042,103.294782318460,103.927027902737,104.556803533993,105.184066658235,105.808775356268,106.430888354220,107.050365033694,107.667165441552,108.281250299330,108.892581012284,109.501119678081,110.106829095126,110.709672770539,111.309614927775,111.906620513896,112.500655206502,113.091685420315,113.679678313426,114.264601793202,114.846424521866,115.425115921748,116.000646180198,116.572986254192,117.142107874608,117.707983550185,118.270586571173,118.829891012672,119.385871737657,119.938504399712,120.487765445448,121.033632116635,121.576082452036,122.115095288947,122.650650264452,123.182727816395,123.711309184063,124.236376408602,124.757912333146,125.275900602683,125.790325663646,126.301172763243,126.808427948523,127.312078065180,127.812110756107,128.308514459691,128.801278407862,129.290392623897,129.775847919974,130.257635894497,130.735748929171,131.210180185857,131.680923603183,132.147973892942,132.611326536254,133.070977779513,133.526924630117,133.979164851980,134.427696960838,134.872520219340,135.313634631940,135.751040939588,136.184740614217,136.614735853038,137.041029572647,137.463625402936,137.882527680821,138.297741443792,138.709272423278,139.117127037838,139.521312386181,139.921836240015,140.318707036728,140.711933871908,141.101526491705,141.487495285030,141.869851275606,142.248606113867,142.623772068712,142.995362019107,143.363389445557,143.727868421433,144.088813604168,144.446240226322,144.800164086515,145.150601540241,145.497569490555,145.841085378648,146.181167174300,146.517833366226,146.851102952309,147.180995429735,147.507530785013,147.830729483911,148.150612461283,148.467201110813,148.780517274661,149.090583233030,149.397421693646,149.701055781157,150.001509026459,150.298805355945,150.592969080687,150.884024885548,151.171997818229,151.456913278262,151.738797005935,152.017675071174,152.293573862364,152.566520075127,152.836540701058,153.103663016413,153.367914570760,153.629323175602,153.887916892958,154.143724023921,154.396773097191,154.647092857583,154.894712254515,155.139660430487,155.381966709535,155.621660585691,155.858771711421,156.093329886071,156.325365044311,156.554907244575,156.781986657520,157.006633554481,157.228878295950,157.448751320066,157.666283131126,157.881504288112,158.094445393255,158.305137080622,158.513610004729,158.719894829204,158.924022215472,159.126022811494,159.325927240545,159.523766090043,159.719569900423,159.913369154077,160.105194264339,160.295075564543,160.483043297130,160.669127602842,160.853358509965,161.035765923664,161.216379615384,161.395229212335,161.572344187059,161.747753847085,161.921487324671,162.093573566644,162.264041324329,162.432919143583,162.600235354928,162.766018063796,162.930295140869,163.093094212551,163.254442651536,163.414367567507,163.572895797950,163.730053899101,163.885868137008,164.040364478739,164.193568583713,164.345505795174,164.496201131806,164.645679279493,164.793964583217,164.941081039122,165.087052286724,165.231901601268,165.375651886264,165.518325666175,165.659945079269,165.800531870652,165.940107385461,166.078692562241,166.216307926502,166.352973584454,166.488709216927,166.623534073489,166.757466966743,166.890526266829,167.022729896124,167.154095324136,167.284639562611,167.414379160845,167.543330201208,167.671508294883,167.798928577817,167.925605706902,168.051553856376,168.176786714448,168.301317480160,168.425158860477,168.548323067616,168.670821816618,168.792666323157,168.913867301603,169.034434963328,169.154379015269,169.273708658746,169.392432588540,169.510558992229,169.628095549799,169.745049433515,169.861427308068,169.977235330998,170.092479153392,170.207163920869,170.321294274846,170.434874354091,170.547907796573,170.660397741601,170.772346832264,170.883757218166,170.994630558481,171.104968025288,171.214770307246,171.324037613568,171.432769678302,171.540965764956,171.648624671432,171.755744735279,171.862323839299,171.968359417457,172.073848461156,172.178787525828,172.283172737881,172.386999801987,172.490264008716,172.592960242537,172.695082990148,172.796626349190,172.897584037310,172.997949401589,173.097715428344,173.196874753294,173.295419672110,173.393342151331,173.490633839668,173.587286079698,173.683289919927,173.778636127266,173.873315199885,173.967317380465,174.060632669862,174.153250841160,174.245161454142,174.336353870157,174.426817267419,174.516540656706,174.605512897482,174.693722714447,174.781158714511,174.867809404190,174.953663207445,175.038708483945,175.122933547778,175.206326686601,175.288876181229,175.370570325686,175.451397447685,175.531345929594,175.610404229823,175.688560904704,175.765804630817,175.842124227791,175.917508681568,175.991947168149,176.065429077815,176.137944039821,176.209481947590,176.280032984367,176.349587649395,176.418136784554,176.485671601518,176.552183709399,176.617665142898,176.682108390961,176.745506425952,176.807852733314,176.869141341774,176.929366854057,176.988524478107,177.046610058853,177.103620110494,177.159551849323,177.214403227064,177.268172964775,177.320860587273,177.372466458110,177.422991815102,177.472438806388,177.520810527073,177.568111056397,177.614345495490,177.659520005673,177.703641847336,177.746719419380,177.788762299231,177.829781283441,177.869788428846,177.908797094337,177.946821983189,177.983879186000,178.019986224199,178.055162094178,178.089427312001,178.122803958718,178.155315726292,178.186987964137,178.217847726239,178.247923818939,178.277246849285,178.305849274037,178.333765449295,178.361031680727,178.387686274464,178.413769588597,178.439324085343,178.464394383814,178.489027313460,178.513271968152,178.537179760895,178.560804479227,178.584202341244,178.607432052297,178.630554862363,178.653634624057,178.676737851338,178.699933778862,178.723294422040,178.746894637742,178.770812185712,178.795127790650,178.819925204986,178.845291272354,178.871315991745,178.898092582372,178.925717549221,178.954290749331,178.983915458758,179.014698440259,179.046750011670,179.080184115053,179.115118386514,179.151674226741,179.189976872321,179.230155467733,179.272343138098,179.316677062660,179.363298549013,179.412353108061,179.463990529730,179.518364959437,179.575634975276,179.635963666016,179.699518709804,179.766472453657,179.837001993729,179.911289256319,179.989521079654,180.071889296475,180.158590817378,180.249827714922,180.345807308564,180.446742250337,180.552850611343,180.664355969039,180.781487495317,180.904480045369,181.033574247391,181.169016593049,181.311059528795];
% magCurveX = [3.91654685730311,3.92068290639734,3.92454735626787,3.92817984882011,3.93161528515215,3.93488442169570,3.93801438034417,3.94102908668770,3.94394964788667,3.94679467964901,3.94958059011725,3.95232182713276,3.95503109425832,3.95771954005548,3.96039692438824,3.96307176492901,3.96575146655085,3.96844243588184,3.97115018295833,3.97387941163011,3.97663410013292,3.97941757304375,3.98223256566581,3.98508128174707,3.98796544531532,3.99088634730892,3.99384488759488,3.99684161288981,3.99987675103482,4.00295024201952,4.00606176610196,4.00921076933011,4.01239648673381,4.01561796342527,4.01887407381848,4.02216353915437,4.02548494349737,4.02883674835114,4.03221730602486,4.03562487186761,4.03905761547583,4.04251363096793,4.04599094641043,4.04948753247139,4.05300131036930,4.05653015917896,4.06007192254950,4.06362441488482,4.06718542703157,4.07075273151565,4.07432408736440,4.07789724454815,4.08146994807183,4.08503994174447,4.08860497165195,4.09216278935622,4.09571115484202,4.09924783923045,4.10277062727711,4.10627731967076,4.10976573514758,4.11323371243436,4.11667911203325,4.12009981785952,4.12349373874273,4.12685880980117,4.13019299369841,4.13349428179028,4.13676069516976,4.13999028561693,4.14318113646041,4.14633136335629,4.14943911499010,4.15250257370704,4.15551995607514,4.15848951338587,4.16140953209622,4.16427833421613,4.16709427764486,4.16985575645948,4.17256120115861,4.17520907886439,4.17779789348513,4.18032618584126,4.18279253375696,4.18519555211945,4.18753389290814,4.18980624519549,4.19201133512124,4.19414792584185,4.19621481745656,4.19821084691154,4.20013488788361,4.20198585064465,4.20376268190797,4.20546436465783,4.20708991796299,4.20863839677554,4.21010889171565,4.21150052884338,4.21281246941820,4.21404390964711,4.21519408042196,4.21626224704688,4.21724770895620,4.21814979942368,4.21896788526357,4.21970136652401,4.22034967617329,4.22091227977955,4.22138867518426,4.22177839216995,4.22208099212268,4.22229606768955,4.22242324243163,XXX4.22246217047269,4.22241253614411,4.22227405362610,4.22204646658578,4.22172954781217,4.22132309884855,4.22082694962227,4.22024095807242,4.21956500977542,4.21879901756892,4.21794292117404,4.21699668681630,4.21596030684535,4.21483379935365,4.21361720779437,4.21231060059857,4.21091407079185,4.20942773561061,4.20785173611815,4.20618623682051,4.20443142528248,4.20258751174362,4.20065472873464,4.19863333069403,4.19652359358523,4.19432581451436,4.19204031134857,4.18966742233514,4.18720750572149,4.18466093937601,4.18202812040993,4.17930946480030,4.17650540701406,4.17361639963335,4.17064291298211,4.16758543475399,4.16444446964174,4.16122053896796,4.15791418031744,4.15452594717111,4.15105640854148,4.14750614860988,4.14387576636539,4.14016587524547,4.13637710277848,4.13251009022802,4.12856549223914,4.12454397648647,4.12044622332440,4.11627292543912,4.11202478750282,4.10770252582986,4.10330686803509,4.09883855269431,4.09429832900679,4.08968695646008,4.08500520449700,4.08025385218479,4.07543368788662,4.07054550893531,4.06559012130941,4.06056833931154,4.05548098524915,4.05032888911755,4.04511288828548,4.03983382718289,4.03449255699131,4.02908993533658,4.02362682598405,4.01810409853624,4.01252262813303,4.00688329515431,4.00118698492519,3.99543458742369,3.98962699699102,3.98376511204441,3.97784983479249,3.97188207095328,3.96586272947475,3.95979272225801,3.95367296388312,3.94750437133747,3.94128786374688,3.93502436210929,3.92871478903112,3.92236006846633,3.91596112545808,3.90951888588315,3.90303427619902,3.89650822319361,3.88994165373781,3.88333549454067,3.87669067190730,3.87000811149953,3.86328873809929,3.85653347537476,3.84974324564920,3.84291896967262,3.83606156639615,3.82917195274919,3.82225104341933,3.81529975063506,3.80831898395125,3.80130965003742,3.79427265246881,3.78720889152021,3.78011926396262,3.77300466286276,3.76586597738527,3.75870409259786,3.75151988927917,3.74431424372950,3.73708802758438,3.72984210763093,3.72257734562708,3.71529459812361,3.70799471628904,3.70067854573735,3.69334692635856,3.68600069215216,3.67864067106338,3.67126768482231,3.66388254878588,3.65648607178274,3.64907905596094,3.64166229663847,3.63423658215680,3.62680269373708,3.61936140533941,3.61191348352486,3.60445968732037,3.59700076808665,3.58953746938879,3.58207052686991,3.57460066812756,3.56712861259310,3.55965507141396,3.55218074733872,3.54470633460519,3.53723251883126,3.52975997690879,3.52228937690023,3.51482137793830,3.50735663012851,3.49989577445446,3.49243944268629,3.48498825729183,3.47754283135070,3.47010376847140,3.46267166271125,3.45524709849918,3.44783065056156,3.44042288385085,3.43302435347715,3.42563560464279,3.41825717257969,3.41088958248968,3.40353334948783,3.39618897854854,3.38885696445469,3.38153779174964,3.37423193469216,3.36693985721427,3.35966201288211,3.35239884485952,3.34515078587479,3.33791825819014,3.33070167357428,3.32350143327773,3.31631792801125,3.30915153792710,3.30200263260319,3.29487157103029,3.28775870160200,3.28066436210788,3.27358887972926,3.26653257103818,3.25949574199913,3.25247868797384,3.24548169372895,3.23850503344655,3.23154897073780,3.22461375865935,3.21769963973280,3.21080684596703,3.20393559888351,3.19708610954451,3.19025857858432,3.18345319624331,3.17667014240500,3.16990958663608,3.16317168822932,3.15645659624947,3.14976444958206,3.14309537698516,3.13644949714414,3.12982691872925,3.12322774045624,3.11665205114998,3.11009992981078,3.10357144568396,3.09706665833219,3.09058561771078,3.08412836424597,3.07769492891615,3.07128533333603,3.06489958984368,3.05853770159069,3.05219966263511,3.04588545803740,3.03959506395935,3.03332844776594,3.02708556813012,3.02086637514052,3.01467081041224,3.00849880720041,3.00235029051684,2.99622517724952,2.99012337628516,2.98404478863465,2.97798930756138,2.97195681871269,2.96594720025405,2.95996032300651,2.95399605058669,2.94805423955004,2.94213473953698,2.93623739342189,2.93036203746517,2.92450850146824,2.91867660893132,2.91286617721447,2.90707701770138,2.90130893596602,2.89556173194261,2.88983520009802,2.88412912960775,2.87844330453417,2.87277750400838,2.86713150241452,2.86150506957727,2.85589797095240,2.85030996781994,2.84474081748061,2.83919027345513,2.83365808568639,2.82814400074468,2.82264776203582,2.81716911001231,2.81170778238733,2.80626351435176,2.80083603879422,2.79542508652389,2.79003038649655,2.78465166604325,2.77928865110222,2.77394106645358,2.76860863595709,2.76329108279280,2.75798812970464,2.75269949924705,2.74742491403455,2.74216409699415,2.73691677162082,2.73168266223603,2.72646149424897,2.72125299442094,2.71605689113262,2.71087291465438,2.70570079741938,2.70054027429980,2.69539108288597,2.69025296376834,2.68512566082272,2.68000892149803,2.67490249710742,2.66980614312216,2.66471961946847,2.65964269082735,2.65457512693739,2.64951670290051,2.64446719949065,2.63942640346548,2.63439410788096,2.62937011240900,2.62435422365790,2.61934625549597,2.61434602937786,2.60935337467413,2.60436812900346,2.59939013856821,2.59441925849246,2.58945535316342,2.58449829657571,2.57954797267836,2.57460427572506,2.56966711062724,2.56473639331019,2.55981205107193,2.55489402294527,2.54998226006280,2.54507672602468,2.54017739726960,2.53528426344850,2.53039732780145,2.52551660753728,2.52064213421638,2.51577395413626,2.51091212872027,2.50605673490910,2.50120786555537,2.49636562982113,2.49153015357823,2.48670157981188,2.48188006902693,2.47706579965734,2.47225896847822,2.46745979102142,2.46266850199355,2.45788535569720,2.45311062645511,2.44834460903734,2.44358761909114,2.43883999357428,2.43410209119073,2.42937429282979,2.42465700200805,2.41995064531406,2.41525567285636,2.41057255871414,2.40590180139116,2.40124392427226,2.39659947608316,2.39196903135317,2.38735319088059,2.38275258220151,2.37816786006113,2.37359970688833,2.36904883327316,2.36451597844716,2.36000191076689,2.35550742819997,2.35103335881472,2.34658056127210,2.34214992532118,2.33774237229714,2.33335885562251,2.32900036131123,2.32466790847576,2.32036254983699,2.31608537223729,2.31183749715651,2.30762008123090,2.30343431677490,2.29928143230582,2.29516269307199,2.29107940158339,2.28703289814508,2.28302456139429,2.27905580883988,2.27512809740487,2.27124292397212,2.26740182593278,2.26360638173784,2.25985821145253,2.25615897731388,2.25251038429075,2.24891418064759,2.24537215851040,2.24188615443611,2.23845804998493,2.23508977229539,2.23178329466243,2.22854063711862,2.22536386701829,2.22225509962436,2.21921649869853,2.21625027709412,2.21335869735196,2.21054407229941,2.20780876565213,2.20515519261879,2.20258582050905,2.20010316934404,2.19770981247024];
% magCurveV = [21.5410077151671,22.1687511430864,22.7961934111071,23.4234451706698,24.0506125409357,24.6777971897729,25.3050964139751,25.9326032187156,26.5604063962388,27.1885906037945,27.8172364408158,28.4464205253453,29.0762155697127,29.7066904554666,30.3379103075632,30.9699365678170,31.6028270676142,32.2366360998939,32.8714144903996,33.5072096682038,34.1440657355100,34.7820235367341,35.4211207268691,36.0613918391361,36.7028683519251,37.3455787550281,37.9895486151685,38.6348006408296,39.2813547463851,39.9292281155354,40.5784352640528,41.2289881018386,41.8808959942947,42.5341658230146,43.1888020457948,43.8448067559723,44.5021797410886,45.1609185408860,45.8210185046372,46.4824728478138,47.1452727080941,47.8094072007159,48.4748634731759,49.1416267592795,49.8096804325446,50.4790060589615,51.1495834491131,51.8213907096584,52.4944042941820,53.1685990534134,53.8439482848187,54.5204237815691,55.1979958808873,55.8766335117773,56.5563042421395,57.2369743252741,57.9186087457769,58.6011712648305,59.2846244648930,59.9689297937890,60.6540476082052,61.3399372165937,62.0265569214870,62.7138640612265,63.4018150511094,64.0903654239562,64.7794698701014,65.4690822768125,66.1591557671387,66.8496427381932,67.5404948988722,68.2316633070144,68.9230984060028,69.6147500608141,70.3065675935158,70.9984998182184,71.6904950754813,72.3825012661796,73.0744658848324,73.7663360523968,74.4580585485307,75.1495798433278,75.8408461285272,76.5318033482022,77.2223972289298,77.9125733094461,78.6022769697877,79.2914534599270,79.9800479278993,80.6680054474292,81.3552710450570,82.0417897267693,82.7275065041364,83.4123664199611,84.0963145734396,84.7792961448405,85.4612564197031,86.1421408125591,86.8218948901801,87.5004643943562,88.1777952642053,88.8538336580201,89.5285259746537,90.2018188744469,90.8736592997024,91.5439944947068,92.2127720253048,92.8799397980289,93.5454460787869,94.2092395111118,94.8712691339758,95.5314843991723,96.1898351882700,96.8462718291398,97.5007451120610,98.1532063054068,98.8036071709152,99.4518999785459,100.098037520928,100.741973127402,101.383660677656,102.023054614968,102.660109959042,103.294782318460,103.927027902737,104.556803533993,105.184066658235,105.808775356268,106.430888354220,107.050365033694,107.667165441552,108.281250299330,108.892581012284,109.501119678081,110.106829095126,110.709672770539,111.309614927775,111.906620513896,112.500655206502,113.091685420315,113.679678313426,114.264601793202,114.846424521866,115.425115921748,116.000646180198,116.572986254192,117.142107874608,117.707983550185,118.270586571173,118.829891012672,119.385871737657,119.938504399712,120.487765445448,121.033632116635,121.576082452036,122.115095288947,122.650650264452,123.182727816395,123.711309184063,124.236376408602,124.757912333146,125.275900602683,125.790325663646,126.301172763243,126.808427948523,127.312078065180,127.812110756107,128.308514459691,128.801278407862,129.290392623897,129.775847919974,130.257635894497,130.735748929171,131.210180185857,131.680923603183,132.147973892942,132.611326536254,133.070977779513,133.526924630117,133.979164851980,134.427696960838,134.872520219340,135.313634631940,135.751040939588,136.184740614217,136.614735853038,137.041029572647,137.463625402936,137.882527680821,138.297741443792,138.709272423278,139.117127037838,139.521312386181,139.921836240015,140.318707036728,140.711933871908,141.101526491705,141.487495285030,141.869851275606,142.248606113867,142.623772068712,142.995362019107,143.363389445557,143.727868421433,144.088813604168,144.446240226322,144.800164086515,145.150601540241,145.497569490555,145.841085378648,146.181167174300,146.517833366226,146.851102952309,147.180995429735,147.507530785013,147.830729483911,148.150612461283,148.467201110813,148.780517274661,149.090583233030,149.397421693646,149.701055781157,150.001509026459,150.298805355945,150.592969080687,150.884024885548,151.171997818229,151.456913278262,151.738797005935,152.017675071174,152.293573862364,152.566520075127,152.836540701058,153.103663016413,153.367914570760,153.629323175602,153.887916892958,154.143724023921,154.396773097191,154.647092857583,154.894712254515,155.139660430487,155.381966709535,155.621660585691,155.858771711421,156.093329886071,156.325365044311,156.554907244575,156.781986657520,157.006633554481,157.228878295950,157.448751320066,157.666283131126,157.881504288112,158.094445393255,158.305137080622,158.513610004729,158.719894829204,158.924022215472,159.126022811494,159.325927240545,159.523766090043,159.719569900423,159.913369154077,160.105194264339,160.295075564543,160.483043297130,160.669127602842,160.853358509965,161.035765923664,161.216379615384,161.395229212335,161.572344187059,161.747753847085,161.921487324671,162.093573566644,162.264041324329,162.432919143583,162.600235354928,162.766018063796,162.930295140869,163.093094212551,163.254442651536,163.414367567507,163.572895797950,163.730053899101,163.885868137008,164.040364478739,164.193568583713,164.345505795174,164.496201131806,164.645679279493,164.793964583217,164.941081039122,165.087052286724,165.231901601268,165.375651886264,165.518325666175,165.659945079269,165.800531870652,165.940107385461,166.078692562241,166.216307926502,166.352973584454,166.488709216927,166.623534073489,166.757466966743,166.890526266829,167.022729896124,167.154095324136,167.284639562611,167.414379160845,167.543330201208,167.671508294883,167.798928577817,167.925605706902,168.051553856376,168.176786714448,168.301317480160,168.425158860477,168.548323067616,168.670821816618,168.792666323157,168.913867301603,169.034434963328,169.154379015269,169.273708658746,169.392432588540,169.510558992229,169.628095549799,169.745049433515,169.861427308068,169.977235330998,170.092479153392,170.207163920869,170.321294274846,170.434874354091,170.547907796573,170.660397741601,170.772346832264,170.883757218166,170.994630558481,171.104968025288,171.214770307246,171.324037613568,171.432769678302,171.540965764956,171.648624671432,171.755744735279,171.862323839299,171.968359417457,172.073848461156,172.178787525828,172.283172737881,172.386999801987,172.490264008716,172.592960242537,172.695082990148,172.796626349190,172.897584037310,172.997949401589,173.097715428344,173.196874753294,173.295419672110,173.393342151331,173.490633839668,173.587286079698,173.683289919927,173.778636127266,173.873315199885,173.967317380465,174.060632669862,174.153250841160,174.245161454142,174.336353870157,174.426817267419,174.516540656706,174.605512897482,174.693722714447,174.781158714511,174.867809404190,174.953663207445,175.038708483945,175.122933547778,175.206326686601,175.288876181229,175.370570325686,175.451397447685,175.531345929594,175.610404229823,175.688560904704,175.765804630817,175.842124227791,175.917508681568,175.991947168149,176.065429077815,176.137944039821,176.209481947590,176.280032984367,176.349587649395,176.418136784554,176.485671601518,176.552183709399,176.617665142898,176.682108390961,176.745506425952,176.807852733314,176.869141341774,176.929366854057,176.988524478107,177.046610058853,177.103620110494,177.159551849323,177.214403227064,177.268172964775,177.320860587273,177.372466458110,177.422991815102,177.472438806388,177.520810527073,177.568111056397,177.614345495490,177.659520005673,177.703641847336,177.746719419380,177.788762299231,177.829781283441,177.869788428846,177.908797094337,177.946821983189,177.983879186000,178.019986224199,178.055162094178,178.089427312001,178.122803958718,178.155315726292,178.186987964137,178.217847726239,178.247923818939,178.277246849285,178.305849274037,178.333765449295,178.361031680727,178.387686274464,178.413769588597,178.439324085343,178.464394383814,178.489027313460,178.513271968152,178.537179760895,178.560804479227,178.584202341244,178.607432052297,178.630554862363,178.653634624057,178.676737851338,178.699933778862,178.723294422040,178.746894637742,178.770812185712,178.795127790650,178.819925204986,178.845291272354,178.871315991745,178.898092582372,178.925717549221,178.954290749331,178.983915458758,179.014698440259,179.046750011670,179.080184115053,179.115118386514,179.151674226741,179.189976872321,179.230155467733,179.272343138098,179.316677062660,179.363298549013,179.412353108061,179.463990529730,179.518364959437,179.575634975276,179.635963666016,179.699518709804,179.766472453657,179.837001993729,179.911289256319,179.989521079654,180.071889296475,180.158590817378,180.249827714922,180.345807308564,180.446742250337,180.552850611343,180.664355969039,180.781487495317,180.904480045369,181.033574247391,181.169016593049,181.311059528795];


Vline_rated = 220;
freq = 60;
f_rated = freq;
poles = 4;
n_mech_init = 1800;
n_mech_max = 1900;
n_mech_min = 1550;

Rs = 0.148/5.3;
Ls = 0.423/(2*pi*f_rated)/5.3;
Rr = 0.144/5.3;
Lr = 0.252/(2*pi*f_rated)/5.3;
Rm = inf;
Rx = 0;
Cx = 100e-05;

R1 = Rs;
Xs = Ls*(2*pi*f_rated);
L1 = Ls;
R2 = Rr;
Xr = Lr*(2*pi*f_rated);
L2 = Lr;
Xm = .4;
Xx = 1/(2*pi*f_rated*Cx);

K_b = 1;
K_w = .005;
