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

Tin_source = 368;
Tin_sink = 278;
m_dot_source = 12;
m_dot_sink = 12;
m_dot_wf_init = 3;
H_wf_init = 2.3e+05;
p_hi = 1e6;
p_low = 1.5e5;
p_atm = 101325;

U_evap = 1500;
A_evap = 9.3;
U_cond = 1400;
A_cond = 18.6;

eff_turbine = 0.78;
eff_pump = 0.78;
eff_inverter = 0.93;
pf_grid = 0.9;

P_spd_speed = [834.591559026165	844.140855611678	853.690152197192	863.239448782706	872.788745368220	882.338041953733	891.887338539247	901.436635124761	910.985931710274	920.535228295788	930.084524881302	939.633821466816	949.183118052329	958.732414637843	968.281711223357	977.831007808871];
P_spd_power = -1/3*[166663.989372846	163937.628261044	155852.617467788	140628.948970298	116654.771019611	83123.6124212207	40806.0621676137	-7481.19463142690	-57239.4174183274	-103560.341265182	-142633.126036591	-172599.754334336	-193461.262910890	-206409.980979167	-213110.460435838	-215209.127126632];

magCurveX = [21.5410077151772,21.5637559851946,21.5850104594817,21.6049891685182,21.6238840683437,21.6418643193326,21.6590790918986,21.6756599767875,21.6917230633813,21.7073707380737,21.7226932456486,21.7377700492335,21.7526710184237,21.7674574703078,21.7821830841376,21.7968947071116,21.8116330660314,21.8264333973516,21.8413260062721,21.8563367639667,21.8714875507320,21.8867966517414,21.9022791111625,21.9179470496093,21.9338099492346,21.9498749101992,21.9661468817719,21.9826288708939,21.9993221306914,22.0162263311072,22.0333397135605,22.0506592313153,22.0681806770356,22.0858987988386,22.1038074060012,22.1218994653485,22.1401671892350,22.1586021159307,22.1771951831361,22.1959367952712,22.2148168851164,22.2338249703230,22.2529502052567,22.2721814285920,22.2915072070305,22.3109158754836,22.3303955740216,22.3499342818659,22.3695198486730,22.3891400233355,22.4087824805036,22.4284348450142,22.4480847143945,22.4677196795940,22.4873273440852,22.5068953414587,22.5264113516305,22.5458631157670,22.5652384500236,22.5845252581887,22.6037115433112,22.6227854183885,22.6417351161824,22.6605489982270,22.6792155630846,22.6977234539060,22.7160614653409,22.7342185498462,22.7521838234333,22.7699465708928,22.7874962505320,22.8048224984593,22.8219151324453,22.8387641553885,22.8553597584131,22.8716923236221,22.8877524265290,22.9035308381886,22.9190185270466,22.9342066605270,22.9490866063722,22.9636499337541,22.9778884141681,22.9917940221269,23.0053589356632,23.0185755366569,23.0314364109947,23.0439343485751,23.0560623431668,23.0678135921302,23.0791814960111,23.0901596580135,23.1007418833599,23.1109221785456,23.1206947504939,23.1300540056181,23.1389945487965,23.1475111822655,23.1555989044362,23.1632529086387,23.1704685818002,23.1772415030592,23.1835674423209,23.1894423587580,23.1948623992592,23.1998238968303,23.2043233689498,23.2083575158822,23.2119232189532,23.2150175387877,23.2176377135136,23.2197811569349,23.2214454566749,23.2226283722927,23.2233278333741,23.2235419376000,23.2232689487928,23.2225072949437,23.2212555662219,23.2195125129671,23.2172770436671,23.2145482229226,23.2113252693984,23.2076075537650,23.2033945966292,23.1986860664574,23.1934817774898,23.1877816876496,23.1815858964452,23.1748946428692,23.1677083032923,23.1600273893553,23.1518525458585,23.1431845486499,23.1340243025129,23.1243728390537,23.1142313145900,23.1036010080406,23.0924833188173,23.0808797647189,23.0687919798291,23.0562217124172,23.0431708228434,23.0296412814683,23.0156351665681,23.0011546622547,22.9862020564017,22.9707797385774,22.9548901979835,22.9385360214017,22.9217198911470,22.9044445830296,22.8867129643238,22.8685279917460,22.8498927094412,22.8308102469782,22.8112838173544,22.7913167150097,22.7709123138501,22.7500740652817,22.7288054962541,22.7071102073153,22.6849918706756,22.6624542282842,22.6395010899152,22.6161363312655,22.5923638920642,22.5681877741930,22.5436120398187,22.5186408095373,22.4932782605304,22.4675286247335,22.4413961870163,22.4148852833764,22.3880002991442,22.3607456672017,22.3331258662134,22.3051454188703,22.2768088901465,22.2481208855701,22.2190860495058,22.1897090634522,22.1599946443512,22.1299475429122,22.0995725419493,22.0688744547316,22.0378581233487,22.0065284170885,21.9748902308302,21.9429484834506,21.9107081162442,21.8781740913586,21.8453513902430,21.8122450121111,21.7788599724190,21.7452013013571,21.7112740423560,21.6770832506078,21.6426339916010,21.6079313396711,21.5729803765647,21.5377861900194,21.5023538723573,21.4666885190945,21.4307952275648,21.3946790955579,21.3583452199737,21.3217986954901,21.2850446132473,21.2480880595460,21.2109341145611,21.1735878510705,21.1360543331994,21.0983386151788,21.0604457401205,21.0223807388063,20.9841486284928,20.9457544117318,20.9072030752058,20.8684995885784,20.8296489033611,20.7906559517944,20.7515256457451,20.7122628756190,20.6728725092882,20.6333593910354,20.5937283405122,20.5539841517141,20.5141315919701,20.4741754009490,20.4341202896799,20.3939709395897,20.3537320015554,20.3134080949721,20.2730038068369,20.2325236908486,20.1919722665227,20.1513540183224,20.1106733948051,20.0699348077852,20.0291426315116,19.9883012018624,19.9474148155540,19.9064877293668,19.8655241593868,19.8245282802621,19.7835042244766,19.7424560816384,19.7013878977846,19.6603036747016,19.6192073692621,19.5781028927769,19.5369941103630,19.4958848403286,19.4547788535720,19.4136798729984,19.3725915729513,19.3315175786608,19.2904614657069,19.2494267594996,19.2084169347747,19.1674354151052,19.1264855724289,19.0855707265928,19.0446941449120,19.0038590417456,18.9630685780887,18.9223258611798,18.8816339441245,18.8409958255355,18.8004144491884,18.7598927036934,18.7194334221833,18.6790393820172,18.6387133045010,18.5984578546232,18.5582756408071,18.5181692146787,18.4781410708518,18.4381936467275,18.3983293223116,18.3585504200460,18.3188592046587,18.2792578830277,18.2397486040621,18.2003334585993,18.1610144793178,18.1217936406668,18.0826728588113,18.0436539915936,18.0047388385112,17.9659291407102,17.9272265809955,17.8886327838565,17.8501493155095,17.8117776839563,17.7735193390582,17.7353756726267,17.6973480185307,17.6594376528190,17.6216457938596,17.5839736024952,17.5464221822141,17.5089925793385,17.4716857832278,17.4345027264987,17.3974442852616,17.3605112793724,17.3237044727016,17.2870245734187,17.2504722342931,17.2140480530111,17.1777525725097,17.1415862813252,17.1055496139596,17.0696429512621,17.0338666208274,16.9982208974096,16.9627060033533,16.9273221090392,16.8920693333485,16.8569477441406,16.8219573587492,16.7870981444935,16.7523700192061,16.7177728517768,16.6833064627130,16.6489706247160,16.6147650632732,16.5806894572677,16.5467434396027,16.5129265978431,16.4792384748727,16.4456785695689,16.4122463374910,16.3789411915880,16.3457625029201,16.3127096013978,16.2797817765362,16.2469782782272,16.2142983175257,16.1817410674539,16.1493056638209,16.1169912060589,16.0847967580758,16.0527213491227,16.0207639746800,15.9889235973579,15.9571991478136,15.9255895256847,15.8940936005397,15.8627102128430,15.8314381749384,15.8002762720465,15.7692232632803,15.7382778826755,15.7074388402387,15.6767048230100,15.6460744961438,15.6155465040037,15.5851194712756,15.5547920040963,15.5245626911976,15.4944301050682,15.4643928031309,15.4344493289353,15.4045982133687,15.3748379758820,15.3451671257316,15.3155841632384,15.2860875810627,15.2566758654952,15.2273474977645,15.1981009553609,15.1689347133759,15.1398472458594,15.1108370271906,15.0819025334682,15.0530422439150,15.0242546422987,14.9955382183698,14.9668914693157,14.9383129012301,14.9098010305996,14.8813543858071,14.8529715086495,14.8246509558734,14.7963913007265,14.7681911345254,14.7400490682397,14.7119637340914,14.6839337871725,14.6559579070772,14.6280347995510,14.6001631981563,14.5723418659534,14.5445695971992,14.5168452190607,14.4891675933460,14.4615356182501,14.4339482301191,14.4064044052284,14.3789031615789,14.3514435607084,14.3240247095197,14.2966457621258,14.2693059217091,14.2420044423996,14.2147406311671,14.1875138497315,14.1603235164885,14.1331691084506,14.1060501632068,14.0789662808964,14.0519171261997,14.0249024303462,13.9979219931365,13.9709756849835,13.9440634489677,13.9171853029087,13.8903413414558,13.8635317381908,13.8367567477501,13.8100167079622,13.7833120420007,13.7566432605554,13.7300109640169,13.7034158446810,13.6768586889662,13.6503403796491,13.6238618981163,13.5974243266311,13.5710288506187,13.5446767609653,13.5183694563355,13.4921084455039,13.4658953497062,13.4397319050022,13.4136199646595,13.3875615015498,13.3615586105650,13.3356135110452,13.3097285492284,13.2839062007110,13.2581490729288,13.2324599076524,13.2068415834984,13.1812971184583,13.1558296724433,13.1304425498443,13.1051392021093,13.0799232303373,13.0547983878867,13.0297685830035,13.0048378814605,12.9800105092190,12.9552908551009,12.9306834734820,12.9061930869975,12.8818245892675,12.8575830476352,12.8334737059250,12.8095019872130,12.7856734966179,12.7619940241045,12.7384695473062,12.7151062343620,12.6919104467714,12.6688887422632,12.6460478776831,12.6233948118972,12.6009367087097,12.5786809397991,12.5566350876699,12.5348069486206,12.5132045357281,12.4918360818479,12.4707100426314,12.4498350995594,12.4292201629903,12.4088743752278,12.3888071136004,12.3690279935629,12.3495468718083,12.3303738493998,12.3115192749185,12.2929937476258,12.2748081206444,12.2569735041537,12.2395012686019,12.2224030479352,12.2056907428433,12.1893765240187,12.1734728354368,12.1579923976478,12.1429482110879,12.1283535594045,12.1142220128009,12.1005674313932,12.0874039685877];
magCurveV = [21.5410077151772,22.1687511430959,22.7961934111159,23.4234451706781,24.0506125409434,24.6777971897801,25.3050964139817,25.9326032187217,26.5604063962444,27.1885906037997,27.8172364408205,28.4464205253496,29.0762155697167,29.7066904554701,30.3379103075664,30.9699365678198,31.6028270676167,32.2366360998961,32.8714144904015,33.5072096682055,34.1440657355114,34.7820235367353,35.4211207268700,36.0613918391368,36.7028683519256,37.3455787550284,37.9895486151686,38.6348006408296,39.2813547463849,39.9292281155350,40.5784352640523,41.2289881018380,41.8808959942940,42.5341658230137,43.1888020457939,43.8448067559713,44.5021797410876,45.1609185408848,45.8210185046360,46.4824728478125,47.1452727080928,47.8094072007146,48.4748634731745,49.1416267592781,49.8096804325431,50.4790060589600,51.1495834491116,51.8213907096569,52.4944042941805,53.1685990534119,53.8439482848172,54.5204237815676,55.1979958808858,55.8766335117758,56.5563042421381,57.2369743252726,57.9186087457755,58.6011712648292,59.2846244648916,59.9689297937877,60.6540476082039,61.3399372165925,62.0265569214858,62.7138640612253,63.4018150511083,64.0903654239551,64.7794698701003,65.4690822768116,66.1591557671378,66.8496427381923,67.5404948988713,68.2316633070136,68.9230984060021,69.6147500608134,70.3065675935151,70.9984998182177,71.6904950754807,72.3825012661791,73.0744658848319,73.7663360523964,74.4580585485304,75.1495798433274,75.8408461285269,76.5318033482019,77.2223972289296,77.9125733094459,78.6022769697876,79.2914534599269,79.9800479278992,80.6680054474292,81.3552710450571,82.0417897267693,82.7275065041365,83.4123664199612,84.0963145734397,84.7792961448407,85.4612564197034,86.1421408125593,86.8218948901805,87.5004643943565,88.1777952642056,88.8538336580205,89.5285259746541,90.2018188744473,90.8736592997029,91.5439944947073,92.2127720253053,92.8799397980294,93.5454460787874,94.2092395111123,94.8712691339763,95.5314843991729,96.1898351882706,96.8462718291404,97.5007451120616,98.1532063054075,98.8036071709159,99.4518999785465,100.098037520929,100.741973127402,101.383660677657,102.023054614969,102.660109959043,103.294782318461,103.927027902738,104.556803533993,105.184066658235,105.808775356269,106.430888354220,107.050365033694,107.667165441553,108.281250299330,108.892581012284,109.501119678081,110.106829095127,110.709672770540,111.309614927775,111.906620513896,112.500655206503,113.091685420316,113.679678313426,114.264601793202,114.846424521867,115.425115921748,116.000646180198,116.572986254193,117.142107874608,117.707983550185,118.270586571173,118.829891012672,119.385871737658,119.938504399712,120.487765445448,121.033632116636,121.576082452037,122.115095288947,122.650650264453,123.182727816395,123.711309184063,124.236376408602,124.757912333146,125.275900602683,125.790325663646,126.301172763243,126.808427948523,127.312078065180,127.812110756107,128.308514459691,128.801278407862,129.290392623897,129.775847919974,130.257635894497,130.735748929171,131.210180185856,131.680923603183,132.147973892942,132.611326536253,133.070977779513,133.526924630117,133.979164851980,134.427696960838,134.872520219340,135.313634631940,135.751040939588,136.184740614216,136.614735853038,137.041029572647,137.463625402936,137.882527680821,138.297741443792,138.709272423277,139.117127037838,139.521312386181,139.921836240015,140.318707036728,140.711933871908,141.101526491705,141.487495285029,141.869851275605,142.248606113867,142.623772068712,142.995362019107,143.363389445556,143.727868421432,144.088813604168,144.446240226322,144.800164086515,145.150601540240,145.497569490554,145.841085378648,146.181167174300,146.517833366225,146.851102952309,147.180995429734,147.507530785012,147.830729483910,148.150612461283,148.467201110813,148.780517274661,149.090583233030,149.397421693646,149.701055781157,150.001509026458,150.298805355944,150.592969080687,150.884024885548,151.171997818229,151.456913278262,151.738797005935,152.017675071174,152.293573862364,152.566520075127,152.836540701058,153.103663016413,153.367914570760,153.629323175602,153.887916892958,154.143724023921,154.396773097191,154.647092857583,154.894712254516,155.139660430487,155.381966709536,155.621660585691,155.858771711421,156.093329886071,156.325365044311,156.554907244576,156.781986657520,157.006633554481,157.228878295950,157.448751320067,157.666283131126,157.881504288112,158.094445393256,158.305137080622,158.513610004730,158.719894829205,158.924022215473,159.126022811495,159.325927240546,159.523766090044,159.719569900424,159.913369154078,160.105194264340,160.295075564544,160.483043297132,160.669127602843,160.853358509966,161.035765923665,161.216379615385,161.395229212337,161.572344187061,161.747753847087,161.921487324673,162.093573566646,162.264041324331,162.432919143585,162.600235354930,162.766018063797,162.930295140871,163.093094212553,163.254442651538,163.414367567508,163.572895797952,163.730053899102,163.885868137010,164.040364478741,164.193568583715,164.345505795176,164.496201131809,164.645679279495,164.793964583219,164.941081039126,165.087052286726,165.231901601270,165.375651886267,165.518325666178,165.659945079272,165.800531870655,165.940107385464,166.078692562244,166.216307926505,166.352973584456,166.488709216930,166.623534073492,166.757466966746,166.890526266833,167.022729896127,167.154095324139,167.284639562614,167.414379160847,167.543330201212,167.671508294887,167.798928577820,167.925605706905,168.051553856379,168.176786714452,168.301317480164,168.425158860480,168.548323067620,168.670821816621,168.792666323161,168.913867301607,169.034434963332,169.154379015273,169.273708658750,169.392432588543,169.510558992233,169.628095549803,169.745049433519,169.861427308073,169.977235331002,170.092479153397,170.207163920873,170.321294274850,170.434874354095,170.547907796578,170.660397741605,170.772346832268,170.883757218171,170.994630558485,171.104968025293,171.214770307251,171.324037613573,171.432769678306,171.540965764960,171.648624671436,171.755744735284,171.862323839303,171.968359417462,172.073848461161,172.178787525833,172.283172737886,172.386999801992,172.490264008722,172.592960242542,172.695082990153,172.796626349195,172.897584037315,172.997949401594,173.097715428350,173.196874753300,173.295419672115,173.393342151337,173.490633839675,173.587286079704,173.683289919934,173.778636127273,173.873315199891,173.967317380471,174.060632669868,174.153250841166,174.245161454148,174.336353870162,174.426817267427,174.516540656712,174.605512897487,174.693722714452,174.781158714517,174.867809404196,174.953663207452,175.038708483952,175.122933547785,175.206326686607,175.288876181236,175.370570325693,175.451397447692,175.531345929599,175.610404229830,175.688560904711,175.765804630824,175.842124227798,175.917508681575,175.991947168156,176.065429077821,176.137944039829,176.209481947597,176.280032984376,176.349587649402,176.418136784562,176.485671601525,176.552183709406,176.617665142906,176.682108390970,176.745506425960,176.807852733322,176.869141341784,176.929366854066,176.988524478114,177.046610058861,177.103620110504,177.159551849332,177.214403227075,177.268172964784,177.320860587283,177.372466458121,177.422991815112,177.472438806400,177.520810527082,177.568111056406,177.614345495499,177.659520005681,177.703641847345,177.746719419388,177.788762299242,177.829781283450,177.869788428856,177.908797094348,177.946821983202,177.983879186012,178.019986224210,178.055162094189,178.089427312010,178.122803958729,178.155315726303,178.186987964148,178.217847726251,178.247923818952,178.277246849295,178.305849274052,178.333765449308,178.361031680742,178.387686274477,178.413769588612,178.439324085356,178.464394383826,178.489027313473,178.513271968163,178.537179760909,178.560804479241,178.584202341258,178.607432052310,178.630554862379,178.653634624071,178.676737851353,178.699933778877,178.723294422055,178.746894637755,178.770812185727,178.795127790663,178.819925205003,178.845291272370,178.871315991763,178.898092582387,178.925717549236,178.954290749348,178.983915458778,179.014698440276,179.046750011685,179.080184115070,179.115118386529,179.151674226759,179.189976872339,179.230155467751,179.272343138116,179.316677062679,179.363298549029,179.412353108079,179.463990529750,179.518364959457,179.575634975294,179.635963666033,179.699518709819,179.766472453675,179.837001993749,179.911289256335,179.989521079670,180.071889296495,180.158590817396,180.249827714940,180.345807308584,180.446742250353,180.552850611359,180.664355969055,180.781487495335,180.904480045386,181.033574247407,181.169016593063,181.311059528815];

Vline = 220;
freq = 60;
f_rated = freq;
pole_pairs = 4;
Rs = 0.018;
R1 = Rs;
Xs = 0.32;
X1 = Xs;
Ls = Xs/(2*pi*f_rated);
L1 = Ls;
Rr = 0.03;
R2 = Rr;
Xr = 0.05;
X2 = Xr;
Lr = Xr/(2*pi*f_rated);
L2 = Lr;
Rm = 5000;
Xm = 18;
Rx = 0.00;
Xx = 12;
Cx = 1/(Xx*2*pi*f_rated);

