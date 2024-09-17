package stockholm_test

import (
	"strings"
	"testing"

	"github.com/evolbioinfo/goalign/io/stockholm"
)

func TestParser_ParseTree(t *testing.T) {
	goodstockohlm := [...]string{
		`# STOCKHOLM 1.0
#=GF ID   Piwi
#=GF AC   PF02171.22
#=GF DE   Piwi domain
#=GF AU   Bateman A;0000-0002-6982-4660
#=GF AU   Hammonds G;
#=GF SE   Bateman A
#=GF GA   28.90 28.90;
#=GF TC   28.90 28.90;
#=GF NC   28.80 28.80;
#=GF BM   hmmbuild HMM.ann SEED.ann
#=GF SM   hmmsearch --cpu 4 -E 1000 -Z 81514348 HMM pfamseq
#=GF TP   Family
#=GF WK   Argonaute
#=GF WK   Piwi
#=GF CL   CL0219
#=GF RC   Redefinition of the Piwi domain
#=GF RN   [1]
#=GF RM   11050429
#=GF RT   Domains in gene silencing and cell differentiation proteins: the
#=GF RT   novel PAZ domain and redefinition of the Piwi domain.
#=GF RA   Cerutti L, Mian N, Bateman A;
#=GF RL   Trends Biochem Sci 2000;25:481-482.
#=GF RN   [2]
#=GF RM   15284453
#=GF RT   Crystal structure of Argonaute and its implications for RISC
#=GF RT   slicer activity.
#=GF RA   Song JJ, Smith SK, Hannon GJ, Joshua-Tor L;
#=GF RL   Science 2004;305:1434-1437.
#=GF DR   INTERPRO; IPR003165;
#=GF DR   SO; 0100021; polypeptide_conserved_region;
#=GF CC   This domain is found in the protein Piwi and its relatives.  The
#=GF CC   function of this domain is the dsRNA guided hydrolysis of ssRNA.
#=GF CC   Determination of the crystal structure of Argonaute reveals that
#=GF CC   PIWI is an RNase H domain, and identifies Argonaute as Slicer,
#=GF CC   the enzyme that cleaves mRNA in the RNAi RISC complex [2].  In
#=GF CC   addition, Mg+2 dependence and production of 3'-OH and 5'
#=GF CC   phosphate products are shared characteristics of RNaseH and
#=GF CC   RISC. The PIWI domain core has a tertiary structure belonging to
#=GF CC   the RNase H family of enzymes.  RNase H fold proteins all have a
#=GF CC   five-stranded mixed beta-sheet surrounded by helices. By analogy
#=GF CC   to RNase H enzymes which cleave single-stranded RNA guided by
#=GF CC   the DNA strand in an RNA/DNA hybrid, the PIWI domain can be
#=GF CC   inferred to cleave single-stranded RNA, for example mRNA, guided
#=GF CC   by double stranded siRNA.
#=GF SQ   15
#=GS AGO1_SCHPO/500-799    AC O74957.1
#=GS AGO6_ARATH/541-851    AC O48771.2
#=GS AGO4_ARATH/577-885    AC Q9ZVD5.2
#=GS TAG76_CAEEL/660-966   AC P34681.2
#=GS O16720_CAEEL/566-867  AC O16720.2
#=GS WAGO4_CAEEL/594-924   AC O62275.1
#=GS YQ53_CAEEL/650-977    AC Q09249.1
#=GS NRDE3_CAEEL/673-1001  AC Q21691.1
#=GS PIWI_ARCFU/110-406    AC O28951.1
#=GS PIWI_ARCFU/110-406    DR PDB; 6XU0 B; 110-406;
#=GS PIWI_ARCFU/110-406    DR PDB; 1YTU B; 110-406;
#=GS PIWI_ARCFU/110-406    DR PDB; 1W9H A; 110-406;
#=GS PIWI_ARCFU/110-406    DR PDB; 6TUO A; 110-406;
#=GS PIWI_ARCFU/110-406    DR PDB; 2W42 B; 110-406;
#=GS PIWI_ARCFU/110-406    DR PDB; 2BGG A; 110-406;
#=GS PIWI_ARCFU/110-406    DR PDB; 2W42 A; 110-406;
#=GS PIWI_ARCFU/110-406    DR PDB; 6XU0 A; 110-406;
#=GS PIWI_ARCFU/110-406    DR PDB; 1YTU A; 110-406;
#=GS PIWI_ARCFU/110-406    DR PDB; 6XUP B; 110-406;
#=GS PIWI_ARCFU/110-406    DR PDB; 6T5T A; 110-406;
#=GS PIWI_ARCFU/110-406    DR PDB; 2BGG B; 110-406;
#=GS PIWI_ARCFU/110-406    DR PDB; 6XUP A; 110-406;
#=GS AGO_METJA/426-699     AC Q58717.1
#=GS AGO_METJA/426-699     DR PDB; 5G5S A; 426-699;
#=GS AGO_METJA/426-699     DR PDB; 5G5T A; 426-699;
#=GS AGO_AQUAE/419-694     AC O67434.1
#=GS AGO_AQUAE/419-694     DR PDB; 2F8S B; 419-694;
#=GS AGO_AQUAE/419-694     DR PDB; 2F8S A; 419-694;
#=GS AGO_AQUAE/419-694     DR PDB; 1YVU A; 419-694;
#=GS AGO_AQUAE/419-694     DR PDB; 2NUB A; 419-694;
#=GS AGO_AQUAE/419-694     DR PDB; 2F8T B; 419-694;
#=GS AGO_AQUAE/419-694     DR PDB; 2F8T A; 419-694;
#=GS PIWL1_HUMAN/555-847   AC Q96J94.1
#=GS AUB_DROME/555-852     AC O76922.1
#=GS PIWI_DROME/538-829    AC Q9VKM1.1
#=GS PIWI_DROME/538-829    DR PDB; 6KR6 A; 538-829;
#=GS AGO10_ARATH/625-946   AC Q9XGW1.1
#=GS AGO10_ARATH/625-946   DR PDB; 7SWQ A; 625-946;
#=GS AGO10_ARATH/625-946   DR PDB; 7SWF A; 625-946;
#=GS AGO10_ARATH/625-946   DR PDB; 7SVA A; 625-946;
AGO1_SCHPO/500-799               YLFFILDK.NSPEP.YGSIKRVCNTMLGVPSQCAISKHILQS......KPQYCANLGMKINVKVGGIN.CSLIPKSNP....LGNVPTL.........ILGGDVYHPGVGA..........TGVSIASIVASVD.LNGCKYTAVSRSQPRHQEVIEG.MKD............IVVYLLQGFRAMTKQ.QPQRIIYFRDGTSEGQFLSVINDELSQIKEACH....SLSPKYN..PKILVCTTQKRHHARFFIKNKSDG......................DRNGNPLPGTII...EKHVTHPYQYDFYLISHPSLQGVSVPVHYTVLHDEIQMPPDQF.QTL......CYNLCYVYARAT....SAVSLVPPVYYAHLVSNLARYQDV
#=GR AGO1_SCHPO/500-799    pAS   ......................................................................................................*.................................................*.............................................*..................................................................................................................................................................................*...........
AGO6_ARATH/541-851               FILCILPERKTSDI.YGPWKKICLTEEGIHTQCICPIKI.........SDQYLTNVLLKINSKLGGIN.SLLGIEYSYNIPLINKIPTL.........ILGMDVSHGPPGR.........ADVPSVAAVVGSKCWPLISRYRAAVRTQSPRLEMIDSLFQPIENTE..KGDNGIMNELFVEFYRTSRARKPKQIIIFRDGVSESQFEQVLKIEVDQIIKAYQ....RLGESDV..PKFTVIVAQKNHHTKLFQAKGPE...........................NVPAGTVV...DTKIVHPTNYDFYMCAHAGKIGTSRPAHYHVLLDEIGFSPDDL.QNL......IHSLSYVNQRST....TATSIVAPVRYAHLAAAQVAQFTK
#=GR AGO6_ARATH/541-851    pAS   ......................................................................................................*.................................................*.............................................*..................................................................................................................................................................................*...........
AGO4_ARATH/577-885               FILCVLPDKKNSDL.YGPWKKKNLTEFGIVTQCMAPTRQPND........QYLTNLLLKINAKLGGLN.SMLSVERTPAFTVISKVPTI.........ILGMDVSHGSPGQ.........SDVPSIAAVVSSREWPLISKYRASVRTQPSKAEMIESLVKKNGT.....EDDGIIKELLVDFYTSSNKRKPEHIIIFRDGVSESQFNQVLNIELDQIIEACK....LLDANWN..PKFLLLVAQKNHHTKFFQPTSPE...........................NVPPGTII...DNKICHPKNNDFYLCAHAGMIGTTRPTHYHVLYDEIGFSADEL.QEL......VHSLSYVYQRST....SAISVVAPICYAHLAAAQLGTFMK
#=GR AGO4_ARATH/577-885    pAS   ......................................................................................................*.................................................*.............................................*..................................................................................................................................................................................*...........
TAG76_CAEEL/660-966              CIIVVLQS.KNSDI.YMTVKEQSDIVHGIMSQCVLMKNVSRP......TPATCANIVLKLNMKMGGIN..SRIVADKITNKYLVDQPTM.........VVGIDVTHPTQAEM.......RMNMPSVAAIVANVD.LLPQSYGANVKVQKKCRESVVY.LLD............AIRERIITFYRHTKQ.KPARIIVYRDGVSEGQFSEVLREEIQSIRTACL....AIAEDFR..PPITYIVVQKRHHARIFCKYQNDM.....................VGKAKNVPPGTTV...DTGIVSPEGFDFYLCSHYGVQGTSRPARYHVLLDECKFTADEI.QSI......TYGMCHTYGRCT....RSVSIPTPVYYADLVATRARCHVK
#=GR TAG76_CAEEL/660-966   pAS   ......................................................................................................*...............................................................................................*..................................................................................................................................................................................*...........
O16720_CAEEL/566-867             LIVVVLPG..KTPI.YAEVKRVGDTVLGIATQCVQAKNAIRT......TPQTLSNLCLKMNVKLGGVN.SILLPNVRPR...IFNEPVI.........FLGCDITHPAAGD.........TRKPSIAAVVGSMD.AHPSRYAATVRVQQHRQEIITD.LTY............MVRELLVQFYRNTRF.KPARIVVYRDGVSEGQLFNVLQYELRAIREACV....MLESGYQ..PGITFIAVQKRHHTRLFAADKADQ.....................VGKAFNIPPGTTV...DVGITHPTEFDFFLCSHAGIQGTSRPSHYHVLWDDNDLTADEL.QQL......TYQMCHTYVRCT....RSVSIPAPAYYAHLVAFRARYHLV
#=GR O16720_CAEEL/566-867  pAS   ......................................................................................................*.................................................*.............................................*..................................................................................................................................................................................*...........
WAGO4_CAEEL/594-924              TFVFIITD.DSITT.LHQRYKMIEKDTKMIVQDMKLSKALSV..INAGKRLTLENVINKTNVKLGGSN..YVFVDAKKQL.....DSHL.........IIGVGISAPPAGTKYAM.ENKGVLNPNVIGYAYNA..QHNQEFSGDFVLNSASQDTLAP.IED............IVMHSLNEYQKFHDGGLPRRVIVYRTGTSEGNHGSIMAYEIPLARAAMR....DFSPDIQ....LVYIVVSKDHSFRFFKPDLASLASRPQATSSTASRHSAMPAAPKAWDLNIAPGILV...DSIVTNPACKQFFLNSHITLQGTAKTPLYTVLADDAKVSMTAL.EDI......TYKLCHLHQIVG....LPTSLPTPLYVANEYAKRGRNLWN
YQ53_CAEEL/650-977               DILVGIAR.EKKPD.VHDILKYFEESIGLQTIQLCQQTVDKMMGG.QGGRQTIDNVMRKFNLKCGGTNFFVEIPNAVRGKAVCSNNETLRKKLLEHVQFIGFEISHGASRTLFDRSRSQMDGEPSVVGVSYSLT...NSTQLGGFTYLQTQKEYKLQKLDE............FFPKCVRSYKEHSKT.LPTRIVIYRVGAGEGNFNRVKE.EVEEMRRTFD....KIQPGYR..PHLVVIIAQRASHARVFPSCISGN.....................RATDQNIPSGTCV...ENVLTSYGYDEFILSSQTPLIGTVRPCKYTILVNDAKWSKNEL.MHL......TYFRAFGHQVSY....QPPSVPDVLYAAENLAKRGRNNYK
NRDE3_CAEEL/673-1001             TIVFGIIA.EKRPD.MHDILKYFEEKLGQQTIQISSETADKFMRD.HGGKQTIDNVIRKLNPKCGGTNFLIDVPESVGHRVVCNNSAEMRAKLYAKTQFIGFEMSHTGARTRFDIQKVMFDGDPTVVGVAYSL..KHSAQLGGFSYFQESRLHKLTN.LQE............KMQICLNAYEQSSSY.LPETVVVYRVGSGEGDYPQIVN.EVNEMKLAAR....KKKHGYN..PKFLVICTQRNSHIRVFPEHINERG....................KSMEQNVKSGTCV...DVPGASHGYEEFILCCQTPLIGTVKPTKYTIIVNDCRWSKNEI.MNV......TYHLAFAHQVSY....APPAIPNVSYAAQNLAKRGHNNYK
PIWI_ARCFU/110-406               GIMLVLPE.YNTPL.YYKLKSYLINS..IPSQFMRYDILSNRNL.....TFYVDNLLVQFVSKLGGK...PWILNVDPEK....................GSDIIIGTGAT........RIDNVNLFCFAMVFK.KDGTMLWNEISPIVTSSEYLTY.LKS............TIKKVVYGFKKSNPDWDVEKLTLHVSG....KRPKMKDGETKILKETVE....ELKKQEMVSRDVKYAILHLNETHPFWVMGDPN........................NRFHPYEGTKVKLSSKRYLLTLLQPYLKRNGLEMVTPIKPLSVEIVSDN..WTSEEYYHNVHEILDEIYYLSKMNWRGF....RSRNLPVTVNYPKLVAGIIANVNR
#=GR PIWI_ARCFU/110-406    SS    EEEEEESS.S-HHH.HHHHHHHHHTT..S-EEEEEHHHHTTSTH.....HHHHHHHHHHHHHHTT-B...--EEE--TTS....................STSEEEEEEEE........ESSSSEEEEEEEEEE.TTS-EEEEEE--SEEGGGHHHH.HHH............HHHHHHHHHHHHSTT---SEEEEEEES....S--S-HHHHHHHHHHHHH....HHHHTTSS-TT-EEEEEEEEE---EEEEE-TT........................TT--TTTTEEEESSSSEEEEEEE-TCEEETTEEEEEEEEEEEEEEEEE-..--TTTHHHHHHHHHHHHHHHTT--TTSS....S--SS-HHHHHHHHHHHHHHHHHH
AGO_METJA/426-699                CFALIIGKEKYKDNDYYEILKKQLFDLKIISQNILWENWRKDDK.....GYMTNNLLIQIMGKLGIK...YFILDSKTPYDY................IMGLDTGLGIFGN............HRVGGCTVVYDSEGKIRRIQPIETPAPGERLHLP................YVIEYLEN..KANIDMENKNILFLRDG.......FIQNSERNDLKEISK....ELN......SNIEVISIRKNNKYKVFTSDYRI........................GSVFGNDGIFLPHKTPFGSNPVKLSTWLRFNCGNEEGLK...IN..........ESI.MQL......LYDLTKMNYSALYGEGRYLRIPAPIHYADKFVKALGKNWK
#=GR AGO_METJA/426-699     SS    EEEEEEE---GGGGHHHHHHHHHHHTTTEEEEEEEHHHHHT-GC.....CHHHHHHHHHHHHHTT-B...--EESS--S-SE................EEEEEE-------............-SEEEEEEEEETTS-EEEEEEEE-SSTT----HH................HHHHHHHH..TS----SS-EEEEEESS.......---CCCHHHHHHHHH....HHT......-EEEEEEEESS-S--BS-S--SE........................EEEETTEEEEE-S--TT----EEEEEEEEEETTEEEEE-...--..........CHH.HHH......HHHHHT----TTSSSSSS-SS-HHHHHHHHHHHHHHTT--
#=GR AGO_METJA/426-699     pAS   ......................................................................................................*...............................................................................................*..................................................................................................................................................................................*...........
#=GR AGO_METJA/426-699     sAS   .......................................................................................................................................................*.............................................................................................................................................................................................................................................
AGO_AQUAE/419-694                LVIVFLEEYPKVDP.YKSFLLYDFVKRELLKKMIPSQVILNRTLKNENLKFVLLNVAEQVLAKTGNIP..YKLKEIEGKVDA................FVGIDISRITRDG..........KTVNAVAFTKIFNSKGELVRYYLTSYPAFGEKLTEK................AIGDVFSLLEKLGF.KKGSKIVVHRDG.......RLYRDEVAAFK........KYGELYG..YSLELLEIIKRNNPRFFSNEKFI..............................KGYFYKLSEDSVILATYNQVY.......EGTHQPIKVRKVYGE..LPVEVL.CSQ......ILSLTLMNYSSF....QPIKLPATVHYSDKITKLMLRGIE
#=GR AGO_AQUAE/419-694     SS    EEEEEE---S---S.SS--HHHHHHHHHHHHTT--EEEEEHHHHHHS-HHHHHHHHHHHHHHHTT-BS..-EES--SS--SE................EEEEEEEEEECSS..........CEEEEEEEEEEEETTS-EEEEEEEEETSSSTCHHHH................HHHHHHHHHHHTT-.-TT-EEEEEESS.......---HHHHHHHH........HHHHHTT..-EEEEEEEE-SS---EEES-S--..............................TTEEEEETTTEEEE--S---S.......SS-B--EEEEEEE-S..S-HHHH.HHH......HHHHGGGG--TT....S--SS-TTTTTHHHHHHHHHTSSS
#=GR AGO_AQUAE/419-694     pAS   ......................................................................................................*...............................................................................................*..................................................................................................................................................................................*...........
#=GR AGO_AQUAE/419-694     sAS   .......................................................................................................................................................*.............................................................................................................................................................................................................................................
PIWL1_HUMAN/555-847              IVVCLLSS.NRKDK.YDAIKKYLCTDCPTPSQCVVARTLGKQQT....VMAIATKIALQMNCKMGGE...LWRVDIPL.......KLVM.........IVGIDCYHDMTAG...........RRSIAGFVASIN.EGMTRWFSRCIFQDRGQELVDG.LKV............CLQAALRAW.NSCNEYMPSRIIVYRDGVGDGQLKTLVNYEVPQFLDCLK....SIGRGYN..PRLTVIVVKKRVNTRFFAQSGGR.........................LQNPLPGTVI...DVEVTRPEWYDFFIVSQAVRSGSVSPTHYNVIYDNSGLKPDHI.QRL......TYKLCHIYYNWP....GVIRVPAPCQYAHKLAFLVGQSIH
#=GR PIWL1_HUMAN/555-847   pAS   ......................................................................................................*.................................................*.............................................*..................................................................................................................................................................................*...........
AUB_DROME/555-852                IVMVVMRS.PNEEK.YSCIKKRTCVDRPVPSQVVTLKVIAPRQQKPTGLMSIATKVVIQMNAKLMGA...PWQVVIPL.......HGLM.........TVGFDVCHSPKNK...........NKSYGAFVATMDQKESFRYFSTVNEHIKGQELSEQ.MSV............NMACALRSYQEQHRS.LPERILFFRDGVGDGQLYQVVNSEVNTLKDRLDEIYKSAGKQEG..CRMTFIIVSKRINSRYFT.............................GHRNPVPGTVV...DDVITLPERYDFFLVSQAVRIGTVSPTSYNVISDNMGLNADKL.QML......SYKMTHMYYNYS....GTIRVPAVCHYAHKLAFLVAESIN
#=GR AUB_DROME/555-852     pAS   ......................................................................................................*.................................................*.............................................*..................................................................................................................................................................................*...........
PIWI_DROME/538-829               LILCLVPN.DNAER.YSSIKKRGYVDRAVPTQVVTLKTTKNRSL.....MSIATKIAIQLNCKLGYT...PWMIELPL.......SGLM.........TIGFDIAKSTRDR...........KRAYGALIASMDLQQNSTYFSTV.TECSAFDVLANTLWP............MIAKALRQYQHEHRK.LPSRIVFYRDGVSSGSLKQLFEFEVKDIIEKLKTEYARVQLSP...PQLAYIVVTRSMNTRFFLN.............................GQNPPPGTIV...DDVITLPERYDFYLVSQQVRQGTVSPTSYNVLYSSMGLSPEKM.QKL......TYKMCHLYYNWS....GTTRVPAVCQYAKKLATLVGTNLH
#=GR PIWI_DROME/538-829    SS    -EEEEESS.--HHH.HHHHHHHHHTTS---EEEEEHHHHT-TTH.....HHHHHHHHHHHHHHTT--...S-EE----.......SSEE.........EEEEEEEE-SS-T...........T-EEEEEEEE-TTTTS--EEEEE.EE---------SHHH............HHHHHHHHHHHHHSS.--SEEEEEEE-----THHHHHHHHHHHHHHHHHHHHHHTT---...-EEEEEEEEE-SS--EEBT.............................TB---TTEEE...-SSSS-SSS-EEEEE-S--SSS----EEEEEEEE-S---HHHH.HHH......HHHHTT--TTSS....B--SS-HHHHHHHHHHHHHHHTT-
#=GR PIWI_DROME/538-829    sAS   ......................................................................................................*...............................................................................................*..............................................................................................................................................................................................
AGO10_ARATH/625-946              LLLAILPD.NNGSL.YGDLKRICETELGLISQCCLTKHVFKI......SKQYLANVSLKINVKMGGRN.TVLVDAISCRIPLVSDIPTI.........IFGADVTHPENGE.........ESSPSIAAVVASQDWPEVTKYAGLVCAQAHRQELIQDLYKTWQDPVRGTVSGGMIRDLLISFRKATGQ.KPLRIIFYRDGVSEGQFYQVLLYELDAIRKACA....SLEPNYQ..PPVTFIVVQKRHHTRLFANNHRDKN...................STDRSGNILPGTVV...DTKICHPTEFDFYLCSHAGIQGTSRPAHYHVLWDENNFTADGI.QSL......TNNLCYTYARCT....RSVSIVPPAYYAHLAAFRARFYLE
#=GR AGO10_ARATH/625-946   SS    EEEEEE-S.S-TSH.HHHHHHHHHHTT---BEEEEHHHHHS-......-HHHHHHHHHHHHHTTSS-S.-EEHHHHTT-STTTTTSTEE.........EEEEEEE--SSSS.........SSS-EEEEEEEE--SST---EEEEEEEESTT-SS-TTSS--EEETTTCEE--TCHHHHHHHHHHHHSS.--SEEEEEEBS--SSGHHHHHHHHHHHHHHHHH....HHSTT--..--EEEEEEE-SSS--EEES-TTSTT...................S-SSSSSB-STEEE...-SSSS-TTS-EEEEE-S--SSS----EEEEEEEESS---HHHH.HHH......HHHHTT--SS-S....B--SS-HHHHHHHHHHHHGGGC--
#=GR AGO10_ARATH/625-946   pAS   ......................................................................................................*.................................................*.............................................*..................................................................................................................................................................................*...........
#=GC SS_cons                     EEEEEESS.SGHHHHHHHHHHHHHHTTCCCEEEEEHHHHCTCTHHHS.HHHHHHHHHHHHHHHTTSBS.-CEEECCTTTSSETTTSTEE.........EEEEEEEEECSSS........ESSSCEEEEEEEEEETTTSEEEEEEEEECSTTCCHHHHSHHHEEETTTCEE--THHHHHHHHHHHHHSS.TTSEEEEEECS--SSCHHHHHHHHHHHHHHHHHHHHHHHCHTTTS.TEEEEEEEEESSS--EEESETTCTT...................S-SCTCCTTTTEEEESTTSCCCECESEEEEEECSCCCSSEEEEEEEEEEEESS-S-HHHHHHHHHHHHHHHHHHTTGGTTTSSSSSSS-SS-HHHHHHHHHHHHHHCTCC
#=GC seq_cons                    hllslLsc.pss-h.YtslKKhspschul.oQslthcslhpp......hhphlsNlhlKlNsKlGGtN..hhllchshth......shh.........hlGhDloHsssus..........tssSluulluShc.ttsscahuslpspspppEhlps.lps............hltchLpsap+pspp.pPpRIllaRDGsuEGph.pllp.Elspl+cshc....plstsap..PplshIllpKpp+sRhFssstsp........................ht.NssPGThV...DstlspPpphDFaLsSpsshpGTs+Ps+YpVLhD-hthos-cl.QpL......oYpLsahahpss....pslulPsPlaYA+hlAttstpph+
//
EOF
`,
	}

	for i, instockholm := range goodstockohlm {
		align, err := stockholm.NewParser(strings.NewReader(instockholm)).Parse()
		if err != nil {
			t.Errorf("Stockholm parser error %d ERROR: %s\n", i, err.Error())
		} else {
			if align.NbSequences() != 15 {
				t.Errorf("Alignment should have 15 sequences but has %d\n", align.NbSequences())
			}
			if align.Length() != 389 {
				t.Errorf("Alignment should be 389 nt long, but is %d\n", align.Length())
			}
		}
	}
}
