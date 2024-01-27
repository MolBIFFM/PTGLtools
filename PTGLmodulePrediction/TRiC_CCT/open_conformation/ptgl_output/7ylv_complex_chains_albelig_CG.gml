graph [
	creator "PTGLgraphComputation"
	version "3.7.1"
	ignoreLigands 1
	descriptionEdgeWeightTypes [
		absoluteWeight "number of residue-residue contacts"
		additiveLengthNormalization "absolute weight / (length chain 1 + chain 2)"
		multiplicativeLengthNormalization "absolute weight / (length chain 1 * chain 2)"
		lucidMultiplicativeLengthNormalization "multiplicative length normalization / smallest multiplicative weight"
		squareRootAdditiveLengthNormalization "absolute weight / sqrt(length chain 1 * chain 2)"
		logAdditiveLengthNormalization "absolute weight / log(length chain 1 + chain 2)"
		minimumLengthNormalization "absolute weight / min(length chain 1, length chain 2)"
		multiplicativeIndividualNormalization "(absolute weight / (length chain 1)) * (absolute weight / (length chain 2))"
		additiveIndividualNormalization "(absolute weight / (length chain 1)) + (absolute weight / (length chain 2))"
		additiveGyrationNormalization "absolute weight / (radius of gyration chain 1 + radius of gyration chain 2)"
		additiveGyrationIndividualNormalization "(absolute weight / (radius of gyration chain 1)) + (absolute weight / (radius of gyration chain 2))"
	]
	minContactsForEdge 1
	factorLucidmultiplicativelengthnormalization 0.00000339156653066484878700623030772
	node [
		id 0
		label "A"
		molName "T-complex protein 1 subunit alpha"
		chainLength 543.0
		chainGyradius 333.90203364319797
		chainType "polypeptide(L)"
	]
	node [
		id 1
		label "a"
		molName "T-complex protein 1 subunit alpha"
		chainLength 543.0
		chainGyradius 329.70086463330307
		chainType "polypeptide(L)"
	]
	node [
		id 2
		label "B"
		molName "T-complex protein 1 subunit beta"
		chainLength 517.0
		chainGyradius 297.41732125502193
		chainType "polypeptide(L)"
	]
	node [
		id 3
		label "b"
		molName "T-complex protein 1 subunit beta"
		chainLength 517.0
		chainGyradius 298.22359826816523
		chainType "polypeptide(L)"
	]
	node [
		id 4
		label "G"
		molName "T-complex protein 1 subunit gamma"
		chainLength 518.0
		chainGyradius 311.108161452
		chainType "polypeptide(L)"
	]
	node [
		id 5
		label "g"
		molName "T-complex protein 1 subunit gamma"
		chainLength 518.0
		chainGyradius 312.84838118254254
		chainType "polypeptide(L)"
	]
	node [
		id 6
		label "D"
		molName "T-complex protein 1 subunit delta"
		chainLength 521.0
		chainGyradius 299.01463471937
		chainType "polypeptide(L)"
	]
	node [
		id 7
		label "d"
		molName "T-complex protein 1 subunit delta"
		chainLength 521.0
		chainGyradius 297.67351215379506
		chainType "polypeptide(L)"
	]
	node [
		id 8
		label "E"
		molName "T-complex protein 1 subunit epsilon"
		chainLength 528.0
		chainGyradius 300.3743379894917
		chainType "polypeptide(L)"
	]
	node [
		id 9
		label "e"
		molName "T-complex protein 1 subunit epsilon"
		chainLength 528.0
		chainGyradius 303.8398774700819
		chainType "polypeptide(L)"
	]
	node [
		id 10
		label "Z"
		molName "T-complex protein 1 subunit zeta"
		chainLength 534.0
		chainGyradius 316.3002903266562
		chainType "polypeptide(L)"
	]
	node [
		id 11
		label "z"
		molName "T-complex protein 1 subunit zeta"
		chainLength 534.0
		chainGyradius 315.07960840845755
		chainType "polypeptide(L)"
	]
	node [
		id 12
		label "H"
		molName "T-complex protein 1 subunit eta"
		chainLength 519.0
		chainGyradius 301.25507964453703
		chainType "polypeptide(L)"
	]
	node [
		id 13
		label "h"
		molName "T-complex protein 1 subunit eta"
		chainLength 519.0
		chainGyradius 300.2718806428935
		chainType "polypeptide(L)"
	]
	node [
		id 14
		label "Q"
		molName "T-complex protein 1 subunit theta"
		chainLength 543.0
		chainGyradius 302.45742121640865
		chainType "polypeptide(L)"
	]
	node [
		id 15
		label "q"
		molName "T-complex protein 1 subunit theta"
		chainLength 543.0
		chainGyradius 303.7138107721431
		chainType "polypeptide(L)"
	]
	node [
		id 16
		label "P"
		molName "Phosducin-like protein 2"
		chainLength 118.0
		chainGyradius 139.50002749367573
		chainType "polypeptide(L)"
	]
	edge [
		source 0
		target 4
		label "40"
		absoluteWeight 40
		additiveLengthNormalization 0.03770028275212064090480678605089538
		multiplicativeLengthNormalization 0.00014221008696146817693779019745871
		lucidMultiplicativeLengthNormalization 41.93050193050193050193050193046400996
		multiplicativeIndividualNormalization 0.0056884034784587270775116078983482301101417123516571030383185079317676
		additiveIndividualNormalization 0.15088490226611773573099539950368680
		additiveGyrationNormalization 0.06201452365275612928527529954310915
		additiveGyrationIndividualNormalization 0.24836826416194084740734061596597009
		numHelixHelixContacts 1
		numHelixStrandContacts 6
		numHelixCoilContacts 5
		numStrandStrandContacts 7
		numStrandCoilContacts 2
		numCoilCoilContacts 19
	]
	edge [
		source 0
		target 6
		label "37"
		absoluteWeight 37
		additiveLengthNormalization 0.03477443609022556390977443609022556
		multiplicativeLengthNormalization 0.00013078687748097404410699073534038
		lucidMultiplicativeLengthNormalization 38.56238003838771593090211132433967456
		multiplicativeIndividualNormalization 0.0048391144667960396319586572075941228169372541118333845876501839853236
		additiveIndividualNormalization 0.13915723763975638292983814240216612
		additiveGyrationNormalization 0.05845951268075065635678752016774498
		additiveGyrationIndividualNormalization 0.23455070823041591878095709646185932
		numHelixHelixContacts 1
		numHelixStrandContacts 2
		numHelixCoilContacts 11
		numStrandStrandContacts 9
		numStrandCoilContacts 5
		numCoilCoilContacts 9
	]
	edge [
		source 0
		target 13
		label "16"
		absoluteWeight 16
		additiveLengthNormalization 0.01506591337099811676082862523540490
		multiplicativeLengthNormalization 0.00005677443163471330686225458364825
		lucidMultiplicativeLengthNormalization 16.73988439306358381502890173408722451
		multiplicativeIndividualNormalization 0.0009083909061554129097960733383720642181273663405685249647821103765919
		additiveIndividualNormalization 0.06029444639606553188771436783444576
		additiveGyrationNormalization 0.02522967223906028714626229648371143
		additiveGyrationIndividualNormalization 0.10120328939593012110878943545082200
		numHelixHelixContacts 5
		numHelixStrandContacts 0
		numHelixCoilContacts 10
		numStrandStrandContacts 0
		numStrandCoilContacts 0
		numCoilCoilContacts 1
	]
	edge [
		source 0
		target 15
		label "2"
		absoluteWeight 2
		additiveLengthNormalization 0.00184162062615101289134438305709024
		multiplicativeLengthNormalization 0.00000678313306132969757401246061543
		lucidMultiplicativeLengthNormalization 1.99999999999999999999999999999705151
		multiplicativeIndividualNormalization 0.0000135662661226593951480249212308673339912972402823139980125420130304
		additiveIndividualNormalization 0.00736648250460405156537753222836096
		additiveGyrationNormalization 0.00313668491383536876419708311645061
		additiveGyrationIndividualNormalization 0.01257492753347269951391564095773470
		numHelixHelixContacts 0
		numHelixStrandContacts 0
		numHelixCoilContacts 2
		numStrandStrandContacts 0
		numStrandCoilContacts 0
		numCoilCoilContacts 0
	]
	edge [
		source 1
		target 16
		label "2"
		absoluteWeight 2
		additiveLengthNormalization 0.00302571860816944024205748865355522
		multiplicativeLengthNormalization 0.00003121390891781377781939632300153
		lucidMultiplicativeLengthNormalization 9.20338983050847457627118644066952043
		multiplicativeIndividualNormalization 0.0000624278178356275556387926460030589886693510628335986515591347504448
		additiveIndividualNormalization 0.02063239379467490713862096950401099
		additiveGyrationNormalization 0.00426256649030143877505951473355143
		additiveGyrationIndividualNormalization 0.02040301954463387232915928596184590
		numHelixHelixContacts 0
		numHelixStrandContacts 0
		numHelixCoilContacts 2
		numStrandStrandContacts 0
		numStrandCoilContacts 0
		numCoilCoilContacts 0
	]
	edge [
		source 1
		target 5
		label "51"
		absoluteWeight 51
		additiveLengthNormalization 0.04806786050895381715362865221489161
		multiplicativeLengthNormalization 0.00018131786087587192559568250175985
		lucidMultiplicativeLengthNormalization 53.46138996138996138996138996134006474
		multiplicativeIndividualNormalization 0.0092472109046694682053798075897523415524008617931269864971522430085966
		additiveIndividualNormalization 0.19237825038930011305701913436720067
		additiveGyrationNormalization 0.07937134831626054422707779133119632
		additiveGyrationIndividualNormalization 0.31770393660620337004882852859384218
		numHelixHelixContacts 0
		numHelixStrandContacts 10
		numHelixCoilContacts 10
		numStrandStrandContacts 9
		numStrandCoilContacts 3
		numCoilCoilContacts 19
	]
	edge [
		source 1
		target 7
		label "44"
		absoluteWeight 44
		additiveLengthNormalization 0.04135338345864661654135338345864662
		multiplicativeLengthNormalization 0.00015553034078818534974885384743181
		lucidMultiplicativeLengthNormalization 45.85796545105566218809980806137890251
		multiplicativeIndividualNormalization 0.0068433349946801553889495692869994305955044661951269516406683562917325
		additiveIndividualNormalization 0.16548428259862921213278049366744078
		additiveGyrationNormalization 0.07013356239592100297690873324598565
		additiveGyrationIndividualNormalization 0.28126725587218272786901532356844086
		numHelixHelixContacts 1
		numHelixStrandContacts 4
		numHelixCoilContacts 13
		numStrandStrandContacts 9
		numStrandCoilContacts 7
		numCoilCoilContacts 10
	]
	edge [
		source 1
		target 12
		label "18"
		absoluteWeight 18
		additiveLengthNormalization 0.01694915254237288135593220338983051
		multiplicativeLengthNormalization 0.00006387123558905247022003640660429
		lucidMultiplicativeLengthNormalization 18.83236994219653179190751445085070750
		multiplicativeIndividualNormalization 0.0011496822406029444639606553188771437779835850924536135151534506435027
		additiveIndividualNormalization 0.06783125219557372337367866381375148
		additiveGyrationNormalization 0.02852814077312779624876645234910903
		additiveGyrationIndividualNormalization 0.11434497311373360638658846346039349
		numHelixHelixContacts 6
		numHelixStrandContacts 0
		numHelixCoilContacts 11
		numStrandStrandContacts 0
		numStrandCoilContacts 0
		numCoilCoilContacts 1
	]
	edge [
		source 1
		target 14
		label "1"
		absoluteWeight 1
		additiveLengthNormalization 0.00092081031307550644567219152854512
		multiplicativeLengthNormalization 0.00000339156653066484878700623030772
		lucidMultiplicativeLengthNormalization 1.00000000000000000000000000000000000
		multiplicativeIndividualNormalization 0.0000033915665306648487870062303077168334978243100705784995031355032576
		additiveIndividualNormalization 0.00368324125230202578268876611418048
		additiveGyrationNormalization 0.00158188229496328767455637711557937
		additiveGyrationIndividualNormalization 0.00633930290335189446306874627740153
		numHelixHelixContacts 0
		numHelixStrandContacts 0
		numHelixCoilContacts 1
		numStrandStrandContacts 0
		numStrandCoilContacts 0
		numCoilCoilContacts 0
	]
	edge [
		source 2
		target 3
		label "22"
		absoluteWeight 22
		additiveLengthNormalization 0.02127659574468085106382978723404255
		multiplicativeLengthNormalization 0.00008230791390592205440553109181448
		lucidMultiplicativeLengthNormalization 24.26840610724721181941643689038494003
		multiplicativeIndividualNormalization 0.0018107741059302851969216840199185154730647351742870076957899502037121
		additiveIndividualNormalization 0.08510638297872340425531914893617022
		additiveGyrationNormalization 0.03693500442785409094569541026401496
		additiveGyrationIndividualNormalization 0.14774028841776954183672674832274481
		numHelixHelixContacts 7
		numHelixStrandContacts 0
		numHelixCoilContacts 11
		numStrandStrandContacts 0
		numStrandCoilContacts 0
		numCoilCoilContacts 4
	]
	edge [
		source 2
		target 6
		label "42"
		absoluteWeight 42
		additiveLengthNormalization 0.04046242774566473988439306358381503
		multiplicativeLengthNormalization 0.00015592689256265847926729210675795
		lucidMultiplicativeLengthNormalization 45.97488834520728995348181038543184613
		multiplicativeIndividualNormalization 0.0065489294876316561292262684838337230182991346057462772454400665288075
		additiveIndividualNormalization 0.16185211448003950147944920681474772
		additiveGyrationNormalization 0.07041876207217054042030150418536748
		additiveGyrationIndividualNormalization 0.28167706856135212053838112866766733
		numHelixHelixContacts 1
		numHelixStrandContacts 0
		numHelixCoilContacts 11
		numStrandStrandContacts 12
		numStrandCoilContacts 10
		numCoilCoilContacts 8
	]
	edge [
		source 2
		target 7
		label "14"
		absoluteWeight 14
		additiveLengthNormalization 0.01348747591522157996146435452793834
		multiplicativeLengthNormalization 0.00005197563085421949308909736891932
		lucidMultiplicativeLengthNormalization 15.32496278173576331782727012847826487
		multiplicativeIndividualNormalization 0.0007276588319590729032473631648704136686999038450829196939377851698675
		additiveIndividualNormalization 0.05395070482667983382648306893824924
		additiveGyrationNormalization 0.02352582028495512879766323330320398
		additiveGyrationIndividualNormalization 0.09410329858059314438273360552737398
		numHelixHelixContacts 0
		numHelixStrandContacts 0
		numHelixCoilContacts 11
		numStrandStrandContacts 0
		numStrandCoilContacts 1
		numCoilCoilContacts 2
	]
	edge [
		source 2
		target 8
		label "79"
		absoluteWeight 79
		additiveLengthNormalization 0.07559808612440191387559808612440191
		multiplicativeLengthNormalization 0.00028940273137565207197702362112420
		lucidMultiplicativeLengthNormalization 85.33010594337963777035343766476952359
		multiplicativeIndividualNormalization 0.0228628157786765136861848660688119105344792216165523709044018521774808
		additiveIndividualNormalization 0.30242585428755641521598968407479046
		additiveGyrationNormalization 0.13215306499899955004468985682218463
		additiveGyrationIndividualNormalization 0.52862519467693439279089287691431441
		numHelixHelixContacts 4
		numHelixStrandContacts 8
		numHelixCoilContacts 16
		numStrandStrandContacts 13
		numStrandCoilContacts 8
		numCoilCoilContacts 30
	]
	edge [
		source 3
		target 6
		label "12"
		absoluteWeight 12
		additiveLengthNormalization 0.01156069364161849710982658959537572
		multiplicativeLengthNormalization 0.00004455054073218813693351203050227
		lucidMultiplicativeLengthNormalization 13.13568238434493998670908868155153482
		multiplicativeIndividualNormalization 0.0005346064887862576432021443660272426614493033409193004080087022056230
		additiveIndividualNormalization 0.04624346128001128613698548766135649
		additiveGyrationNormalization 0.02009248460195355276140181833578919
		additiveGyrationIndividualNormalization 0.08037007939907248154634145992742654
		numHelixHelixContacts 0
		numHelixStrandContacts 0
		numHelixCoilContacts 8
		numStrandStrandContacts 0
		numStrandCoilContacts 1
		numCoilCoilContacts 3
	]
	edge [
		source 3
		target 7
		label "42"
		absoluteWeight 42
		additiveLengthNormalization 0.04046242774566473988439306358381503
		multiplicativeLengthNormalization 0.00015592689256265847926729210675795
		lucidMultiplicativeLengthNormalization 45.97488834520728995348181038543184613
		multiplicativeIndividualNormalization 0.0065489294876316561292262684838337230182991346057462772454400665288075
		additiveIndividualNormalization 0.16185211448003950147944920681474772
		additiveGyrationNormalization 0.07048196620765522608298536735982480
		additiveGyrationIndividualNormalization 0.28192810507718099802201479766044753
		numHelixHelixContacts 3
		numHelixStrandContacts 3
		numHelixCoilContacts 10
		numStrandStrandContacts 10
		numStrandCoilContacts 9
		numCoilCoilContacts 7
	]
	edge [
		source 3
		target 9
		label "61"
		absoluteWeight 61
		additiveLengthNormalization 0.05837320574162679425837320574162679
		multiplicativeLengthNormalization 0.00022346286853056678975441064415919
		lucidMultiplicativeLengthNormalization 65.88780332336908739229822401963145465
		multiplicativeIndividualNormalization 0.0136312349803645741750190492937108025400035167926850712150518726921048
		additiveIndividualNormalization 0.23351869761444229529335912314635719
		additiveGyrationNormalization 0.10131822051686844894006358878097781
		additiveGyrationIndividualNormalization 0.40530815147423282705940879558118401
		numHelixHelixContacts 0
		numHelixStrandContacts 10
		numHelixCoilContacts 17
		numStrandStrandContacts 11
		numStrandCoilContacts 4
		numCoilCoilContacts 19
	]
	edge [
		source 4
		target 10
		label "48"
		absoluteWeight 48
		additiveLengthNormalization 0.04562737642585551330798479087452471
		multiplicativeLengthNormalization 0.00017352826341590386534206758925860
		lucidMultiplicativeLengthNormalization 51.16463493991583879224328662526114930
		multiplicativeIndividualNormalization 0.0083293566439633855364192442844128237247841742223764695674808034358596
		additiveIndividualNormalization 0.18255173311353086633985510390004772
		additiveGyrationNormalization 0.07650518551977356607661966315836320
		additiveGyrationIndividualNormalization 0.30604170109103054725458569042519593
		numHelixHelixContacts 2
		numHelixStrandContacts 8
		numHelixCoilContacts 8
		numStrandStrandContacts 11
		numStrandCoilContacts 7
		numCoilCoilContacts 12
	]
	edge [
		source 4
		target 11
		label "2"
		absoluteWeight 2
		additiveLengthNormalization 0.00190114068441064638783269961977186
		multiplicativeLengthNormalization 0.00000723034430899599438925281621911
		lucidMultiplicativeLengthNormalization 2.13185978916315994967680360938637264
		multiplicativeIndividualNormalization 0.0000144606886179919887785056324382167143869391060402296357352537127818
		additiveIndividualNormalization 0.00760632221306378609749396266250199
		additiveGyrationNormalization 0.00319393015364335339664693540594690
		additiveGyrationIndividualNormalization 0.01277623452951666213701530120273570
		numHelixHelixContacts 0
		numHelixStrandContacts 0
		numHelixCoilContacts 2
		numStrandStrandContacts 0
		numStrandCoilContacts 0
		numCoilCoilContacts 0
	]
	edge [
		source 4
		target 15
		label "14"
		absoluteWeight 14
		additiveLengthNormalization 0.01319509896324222431668237511781338
		multiplicativeLengthNormalization 0.00004977353043651386192822656911055
		lucidMultiplicativeLengthNormalization 14.67567567567567567567567567566284576
		multiplicativeIndividualNormalization 0.0006968294261111940669951719675476581847593449803394554775770245383505
		additiveIndividualNormalization 0.05280971579314120750584838982629038
		additiveGyrationNormalization 0.02277081924927705411271048122579022
		additiveGyrationIndividualNormalization 0.09109645357996265679451882387849692
		numHelixHelixContacts 1
		numHelixStrandContacts 0
		numHelixCoilContacts 11
		numStrandStrandContacts 0
		numStrandCoilContacts 0
		numCoilCoilContacts 2
	]
	edge [
		source 5
		target 10
		label "3"
		absoluteWeight 3
		additiveLengthNormalization 0.00285171102661596958174904942965779
		multiplicativeLengthNormalization 0.00001084551646349399158387922432866
		lucidMultiplicativeLengthNormalization 3.19778968374473992451520541407808471
		multiplicativeIndividualNormalization 0.0000325365493904819747516376729859875784130840310615591514467918962301
		additiveIndividualNormalization 0.01140948331959567914624094399375298
		additiveGyrationNormalization 0.00476834830280037745217284156009793
		additiveGyrationIndividualNormalization 0.01907396739797062262203458481608578
		numHelixHelixContacts 1
		numHelixStrandContacts 0
		numHelixCoilContacts 2
		numStrandStrandContacts 0
		numStrandCoilContacts 0
		numCoilCoilContacts 0
	]
	edge [
		source 5
		target 11
		label "46"
		absoluteWeight 46
		additiveLengthNormalization 0.04372623574144486692015209125475285
		multiplicativeLengthNormalization 0.00016629791910690787095281477303949
		lucidMultiplicativeLengthNormalization 49.03277515075267884256648301587477667
		multiplicativeIndividualNormalization 0.0076497042789177620638294795598166383053519008575188350469249345653840
		additiveIndividualNormalization 0.17494541090046708024236114123754573
		additiveGyrationNormalization 0.07325680772720774521822999536230102
		additiveGyrationIndividualNormalization 0.29303093073552850157023760818542368
		numHelixHelixContacts 2
		numHelixStrandContacts 8
		numHelixCoilContacts 10
		numStrandStrandContacts 11
		numStrandCoilContacts 4
		numCoilCoilContacts 11
	]
	edge [
		source 5
		target 14
		label "20"
		absoluteWeight 20
		additiveLengthNormalization 0.01885014137606032045240339302544769
		multiplicativeLengthNormalization 0.00007110504348073408846889509872935
		lucidMultiplicativeLengthNormalization 20.96525096525096525096525096523053074
		multiplicativeIndividualNormalization 0.0014221008696146817693779019745870575275354280879142757595796269829419
		additiveIndividualNormalization 0.07544245113305886786549769975184340
		additiveGyrationNormalization 0.03250416284394540058872092049621708
		additiveGyrationIndividualNormalization 0.13005374097762602099399523011044187
		numHelixHelixContacts 2
		numHelixStrandContacts 0
		numHelixCoilContacts 12
		numStrandStrandContacts 0
		numStrandCoilContacts 0
		numCoilCoilContacts 6
	]
	edge [
		source 6
		target 9
		label "25"
		absoluteWeight 25
		additiveLengthNormalization 0.02383222116301239275500476644423260
		multiplicativeLengthNormalization 0.00009088000930611295294596638166696
		lucidMultiplicativeLengthNormalization 26.79587986389809806316524166809645419
		multiplicativeIndividualNormalization 0.0022720002326528238236491595416739370992191589600418775082882568487175
		additiveIndividualNormalization 0.09533312976211248764031873436863840
		additiveGyrationNormalization 0.04146937527133171404883371138008794
		additiveGyrationIndividualNormalization 0.16588812852405795888319994703452176
		numHelixHelixContacts 9
		numHelixStrandContacts 0
		numHelixCoilContacts 13
		numStrandStrandContacts 0
		numStrandCoilContacts 0
		numCoilCoilContacts 3
	]
	edge [
		source 7
		target 8
		label "22"
		absoluteWeight 22
		additiveLengthNormalization 0.02097235462345090562440419447092469
		multiplicativeLengthNormalization 0.00007997440818937939859245041586692
		lucidMultiplicativeLengthNormalization 23.58037428023032629558541266792346441
		multiplicativeIndividualNormalization 0.0017594369801663467690339091490722968074216250799744081893793985924504
		additiveIndividualNormalization 0.08389315419065898912348048624440179
		additiveGyrationNormalization 0.03678635412656195069506948005184493
		additiveGyrationIndividualNormalization 0.14714841758151555566208179194328423
		numHelixHelixContacts 10
		numHelixStrandContacts 0
		numHelixCoilContacts 11
		numStrandStrandContacts 0
		numStrandCoilContacts 0
		numCoilCoilContacts 1
	]
	edge [
		source 8
		target 12
		label "76"
		absoluteWeight 76
		additiveLengthNormalization 0.07258834765998089780324737344794651
		multiplicativeLengthNormalization 0.00027733987271559525894786010392947
		lucidMultiplicativeLengthNormalization 81.77338413032054650551760378342390079
		multiplicativeIndividualNormalization 0.0210778303263852396800373678986395747099608804811116949845273544695510
		additiveIndividualNormalization 0.29037484673322823611840952881415309
		additiveGyrationNormalization 0.12632361013674835300811481997841419
		additiveGyrationIndividualNormalization 0.50529552343712849535782151758428581
		numHelixHelixContacts 8
		numHelixStrandContacts 7
		numHelixCoilContacts 22
		numStrandStrandContacts 11
		numStrandCoilContacts 3
		numCoilCoilContacts 25
	]
	edge [
		source 9
		target 13
		label "81"
		absoluteWeight 81
		additiveLengthNormalization 0.07736389684813753581661891117478510
		multiplicativeLengthNormalization 0.00029558591697320021019442984760904
		lucidMultiplicativeLengthNormalization 87.15321203363110877561744613759640946
		multiplicativeIndividualNormalization 0.0239424592748292170257488176563321066759721492380451918024172359432475
		additiveIndividualNormalization 0.30947845507094062007356805044666316
		additiveGyrationNormalization 0.13408115123104776303197178716782694
		additiveGyrationIndividualNormalization 0.53634331421014971621217853959378263
		numHelixHelixContacts 5
		numHelixStrandContacts 7
		numHelixCoilContacts 28
		numStrandStrandContacts 10
		numStrandCoilContacts 2
		numCoilCoilContacts 29
	]
	edge [
		source 10
		target 11
		label "13"
		absoluteWeight 13
		additiveLengthNormalization 0.01217228464419475655430711610486891
		multiplicativeLengthNormalization 0.00004558908106439983728204912398827
		lucidMultiplicativeLengthNormalization 13.44189496275722762277490215880486274
		multiplicativeIndividualNormalization 0.0005926580538371978846666386118475502241580047412644306975830773331089
		additiveIndividualNormalization 0.04868913857677902621722846441947566
		additiveGyrationNormalization 0.02058982242868957861976767634579136
		additiveGyrationIndividualNormalization 0.08235959756375076246734823569519409
		numHelixHelixContacts 2
		numHelixStrandContacts 0
		numHelixCoilContacts 9
		numStrandStrandContacts 0
		numStrandCoilContacts 0
		numCoilCoilContacts 2
	]
	edge [
		source 10
		target 14
		label "45"
		absoluteWeight 45
		additiveLengthNormalization 0.04178272980501392757660167130919220
		multiplicativeLengthNormalization 0.00015519274939474827736048171829412
		lucidMultiplicativeLengthNormalization 45.75842696629213483146067415726023670
		multiplicativeIndividualNormalization 0.0069836737227636724812216773232354585464647091687876342417282264572599
		additiveIndividualNormalization 0.16714259109814389471723881060276864
		additiveGyrationNormalization 0.07272636633130357348749965648435696
		additiveGyrationIndividualNormalization 0.29105113873256808718435188831802415
		numHelixHelixContacts 4
		numHelixStrandContacts 4
		numHelixCoilContacts 8
		numStrandStrandContacts 11
		numStrandCoilContacts 5
		numCoilCoilContacts 13
	]
	edge [
		source 11
		target 15
		label "54"
		absoluteWeight 54
		additiveLengthNormalization 0.05013927576601671309192200557103064
		multiplicativeLengthNormalization 0.00018623129927369793283257806195295
		lucidMultiplicativeLengthNormalization 54.91011235955056179775280898871405313
		multiplicativeIndividualNormalization 0.0100564901607796883729592153454590605158606989881432739462412316096592
		additiveIndividualNormalization 0.20057110931777267366068657272332237
		additiveGyrationNormalization 0.08726660356457280725046921360073570
		additiveGyrationIndividualNormalization 0.34918421908209262057522560677864124
		numHelixHelixContacts 3
		numHelixStrandContacts 7
		numHelixCoilContacts 12
		numStrandStrandContacts 13
		numStrandCoilContacts 3
		numCoilCoilContacts 16
	]
	edge [
		source 12
		target 14
		label "68"
		absoluteWeight 68
		additiveLengthNormalization 0.06403013182674199623352165725047081
		multiplicativeLengthNormalization 0.00024129133444753155416458198050508
		lucidMultiplicativeLengthNormalization 71.14450867052023121387283236987586403
		multiplicativeIndividualNormalization 0.0164078107424321456831915746743454087021010088106821093120713086861332
		additiveIndividualNormalization 0.25625139718327851052278606329639447
		additiveGyrationNormalization 0.11263639547471053176019449395292917
		additiveGyrationIndividualNormalization 0.45054736894531760210825771732297708
		numHelixHelixContacts 3
		numHelixStrandContacts 9
		numHelixCoilContacts 19
		numStrandStrandContacts 11
		numStrandCoilContacts 7
		numCoilCoilContacts 19
	]
	edge [
		source 13
		target 15
		label "66"
		absoluteWeight 66
		additiveLengthNormalization 0.06214689265536723163841807909604520
		multiplicativeLengthNormalization 0.00023419453049319239080680015754905
		lucidMultiplicativeLengthNormalization 69.05202312138728323699421965311532953
		multiplicativeIndividualNormalization 0.0154568390125506977932488103982371540031296905438635710407817839236100
		additiveIndividualNormalization 0.24871459138377031903682176731708875
		additiveGyrationNormalization 0.10927411185085052643971211297827216
		additiveGyrationIndividualNormalization 0.43711064262511901909723314603227558
		numHelixHelixContacts 5
		numHelixStrandContacts 9
		numHelixCoilContacts 16
		numStrandStrandContacts 12
		numStrandCoilContacts 12
		numCoilCoilContacts 12
	]
	edge [
		source 14
		target 15
		label "4"
		absoluteWeight 4
		additiveLengthNormalization 0.00368324125230202578268876611418048
		multiplicativeLengthNormalization 0.00001356626612265939514802492123087
		lucidMultiplicativeLengthNormalization 3.99999999999999999999999999999705151
		multiplicativeIndividualNormalization 0.0000542650644906375805920996849234693359651889611292559920501680521216
		additiveIndividualNormalization 0.01473296500920810313075506445672192
		additiveGyrationNormalization 0.00659879550350476655517528493288062
		additiveGyrationIndividualNormalization 0.02639529540656823392773229412161272
		numHelixHelixContacts 0
		numHelixStrandContacts 0
		numHelixCoilContacts 0
		numStrandStrandContacts 0
		numStrandCoilContacts 0
		numCoilCoilContacts 4
	]
]
