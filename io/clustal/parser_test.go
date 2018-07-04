package clustal

import (
	"fmt"
	"strings"
	"testing"
)

// Sample file from t-coffee:
// https://github.com/cbcrg/tcoffee/blob/master/examples/sample_aln2.aln
var clustalstring1 string = `CLUSTAL W (1.83) multiple sequence alignment

1cms            --GEVASVPLTNYLDSQYFGKIYLGTPPQEFTVLFDTGSSDFWVPSIYCK
4pep            ----IGDEPLENYLDTEYFGTIGIGTPAQDFTVIFDTGSSNLWVPSVYCS
4ape            -STGSATTTPIDSLDDAYITPVQIGTPAQTLNLDFDTGSSDLWVFSSETT
3app            AASGVATNTPTAN-DEEYITPVTIG--GTTLNLNFDTGSADLWVFSTELP
2apr            AGVGTVPMTDYGN-DIEYYGQVTIGTPGKKFNLDFDTGSSDLWIASTLCT
1cms_1          -----------------YTGSLHWVPVTVQQYWQFTVDS---VTISGVVV
                                 *   :            * ..*      *    

1cms            SNACKNHQRFDPR-KSSTFQNLGKPLSIHYGTG-SMQGILGYDTVTVSNI
4pep            SLACSDHNQFNPD-DSSTFEATSQELSITYGTG-SMTGILGYDTVQVGGI
4ape            ASEVDGQTIYTPSKSTTAKLLSGATWSISYGDGSSSSGDVYTDTVSVGGL
3app            ASQQSGHSVYNP--SATGKELSGYTWSISYGDGSSASGNVFTDSVTVGGV
2apr            NCG-SGQTKYDPN-QSSTYQADGRTWSISYGDGSSASGILAKDNVNLGGL
1cms_1          ACEGGCQAILDTG-TSKLVGPSSDILNIQQAIG-ATQNQYGEFDIDCDNL
                      :    .   :.     .   .*  . * :  .      :  ..:

1cms            VDIQQTVGLSTQEPGDVFTYAEFDGILGMAYPSLASEY---SIPVFDNMM
4pep            SDTNQIFGLSETEPGSFLYYAPFDGILGLAYPSISASG---ATPVFDNLW
4ape            TVTGQAVESAKKVSSSFTEDSTIDGLLGLAFSTLNTVSPTQQKTFFDNAK
3app            TAHGQAVQAAQQISAQFQQDTNNDGLLGLAFSSINTVQPQSQTTFFDTVK
2apr            LIKGQTIELAKREAASFASGPN-DGLLGLGFDTITTVR--GVKTPMDNLI
1cms_1          SYMPTVVFEINGKMYPLTPSAYTSQDQGFCTSGFQSEN------HSQKWI
                      .         .   .  .   *:    : :          :.  

1cms            NRHLVAQDLFSVYMDRNGQE--SMLTLGAIDPSY
4pep            DQGLVSQDLFSVYLSSNDDSG-SVVLLGGIDSSY
4ape            AS--LDSPVFTADLGYHAPG---TYNFGFIDTTA
3app            SS--LAQPLFAVALKHQQPG---VYDFGFIDSSK
2apr            SQGLISRPIFGVYLGKAKNGGGGEYIFGGYDSTK
1cms_1          LGDVFIREYYSVFDRANNLVG----LAKAI----
                    .    : .                      



`

var clustalstring2 string = `CLUSTALW

1cms     --GEVASVPLTNYLDSQYFGKIYLGTPPQE 30
4pep     ----IGDEPLENYLDTEYFGTIGIGTPAQD 30
4ape     -STGSATTTPIDSLDDAYITPVQIGTPAQT 30
3app     AASGVATNTPTAN-DEEYITPVTIG--GTT 30
2apr     AGVGTVPMTDYGN-DIEYYGQVTIGTPGKK 30
1cms_1   -----------------YTGSLHWVPVTVQ 30
                          *   :        


`

// New sequence name in 2nd block
var clustalstring3 string = `CLUSTAL W (1.83) multiple sequence alignment

1cms            --GEVASVPLTNYLDSQYFGKIYLGTPPQEFTVLFDTGSSDFWVPSIYCK
4pep            ----IGDEPLENYLDTEYFGTIGIGTPAQDFTVIFDTGSSNLWVPSVYCS
4ape            -STGSATTTPIDSLDDAYITPVQIGTPAQTLNLDFDTGSSDLWVFSSETT
3app            AASGVATNTPTAN-DEEYITPVTIG--GTTLNLNFDTGSADLWVFSTELP
2apr            AGVGTVPMTDYGN-DIEYYGQVTIGTPGKKFNLDFDTGSSDLWIASTLCT
1cms_1          -----------------YTGSLHWVPVTVQQYWQFTVDS---VTISGVVV
                                 *   :            * ..*      *    

1cms            SNACKNHQRFDPR-KSSTFQNLGKPLSIHYGTG-SMQGILGYDTVTVSNI
4pep            SLACSDHNQFNPD-DSSTFEATSQELSITYGTG-SMTGILGYDTVQVGGI
4ape            ASEVDGQTIYTPSKSTTAKLLSGATWSISYGDGSSSSGDVYTDTVSVGGL
33app            ASQQSGHSVYNP--SATGKELSGYTWSISYGDGSSASGNVFTDSVTVGGV
2apr            NCG-SGQTKYDPN-QSSTYQADGRTWSISYGDGSSASGILAKDNVNLGGL
1cms_1          ACEGGCQAILDTG-TSKLVGPSSDILNIQQAIG-ATQNQYGEFDIDCDNL
                      :    .   :.     .   .*  . * :  .      :  ..:

1cms            VDIQQTVGLSTQEPGDVFTYAEFDGILGMAYPSLASEY---SIPVFDNMM
4pep            SDTNQIFGLSETEPGSFLYYAPFDGILGLAYPSISASG---ATPVFDNLW
4ape            TVTGQAVESAKKVSSSFTEDSTIDGLLGLAFSTLNTVSPTQQKTFFDNAK
3app            TAHGQAVQAAQQISAQFQQDTNNDGLLGLAFSSINTVQPQSQTTFFDTVK
2apr            LIKGQTIELAKREAASFASGPN-DGLLGLGFDTITTVR--GVKTPMDNLI
1cms_1          SYMPTVVFEINGKMYPLTPSAYTSQDQGFCTSGFQSEN------HSQKWI
                      .         .   .  .   *:    : :          :.  

1cms            NRHLVAQDLFSVYMDRNGQE--SMLTLGAIDPSY
4pep            DQGLVSQDLFSVYLSSNDDSG-SVVLLGGIDSSY
4ape            AS--LDSPVFTADLGYHAPG---TYNFGFIDTTA
3app            SS--LAQPLFAVALKHQQPG---VYDFGFIDSSK
2apr            SQGLISRPIFGVYLGKAKNGGGGEYIFGGYDSTK
1cms_1          LGDVFIREYYSVFDRANNLVG----LAKAI----
                    .    : .                      



`

// Lack a sequence in 2nd block
var clustalstring4 string = `CLUSTAL W (1.83) multiple sequence alignment

1cms            --GEVASVPLTNYLDSQYFGKIYLGTPPQEFTVLFDTGSSDFWVPSIYCK
4pep            ----IGDEPLENYLDTEYFGTIGIGTPAQDFTVIFDTGSSNLWVPSVYCS
4ape            -STGSATTTPIDSLDDAYITPVQIGTPAQTLNLDFDTGSSDLWVFSSETT
3app            AASGVATNTPTAN-DEEYITPVTIG--GTTLNLNFDTGSADLWVFSTELP
2apr            AGVGTVPMTDYGN-DIEYYGQVTIGTPGKKFNLDFDTGSSDLWIASTLCT
1cms_1          -----------------YTGSLHWVPVTVQQYWQFTVDS---VTISGVVV
                                 *   :            * ..*      *    

1cms            SNACKNHQRFDPR-KSSTFQNLGKPLSIHYGTG-SMQGILGYDTVTVSNI
4pep            SLACSDHNQFNPD-DSSTFEATSQELSITYGTG-SMTGILGYDTVQVGGI
4ape            ASEVDGQTIYTPSKSTTAKLLSGATWSISYGDGSSSSGDVYTDTVSVGGL
3app            ASQQSGHSVYNP--SATGKELSGYTWSISYGDGSSASGNVFTDSVTVGGV
2apr            NCG-SGQTKYDPN-QSSTYQADGRTWSISYGDGSSASGILAKDNVNLGGL
                      :    .   :.     .   .*  . * :  .      :  ..:

1cms            VDIQQTVGLSTQEPGDVFTYAEFDGILGMAYPSLASEY---SIPVFDNMM
4pep            SDTNQIFGLSETEPGSFLYYAPFDGILGLAYPSISASG---ATPVFDNLW
4ape            TVTGQAVESAKKVSSSFTEDSTIDGLLGLAFSTLNTVSPTQQKTFFDNAK
3app            TAHGQAVQAAQQISAQFQQDTNNDGLLGLAFSSINTVQPQSQTTFFDTVK
2apr            LIKGQTIELAKREAASFASGPN-DGLLGLGFDTITTVR--GVKTPMDNLI
1cms_1          SYMPTVVFEINGKMYPLTPSAYTSQDQGFCTSGFQSEN------HSQKWI
                      .         .   .  .   *:    : :          :.  

1cms            NRHLVAQDLFSVYMDRNGQE--SMLTLGAIDPSY
4pep            DQGLVSQDLFSVYLSSNDDSG-SVVLLGGIDSSY
4ape            AS--LDSPVFTADLGYHAPG---TYNFGFIDTTA
3app            SS--LAQPLFAVALKHQQPG---VYDFGFIDSSK
2apr            SQGLISRPIFGVYLGKAKNGGGGEYIFGGYDSTK
1cms_1          LGDVFIREYYSVFDRANNLVG----LAKAI----
                    .    : .                      



`

func TestParse(t *testing.T) {
	align, err := NewParser(strings.NewReader(clustalstring1)).Parse()

	if err != nil {
		t.Error(err)
	}
	if align.Length() != 184 {
		t.Error(fmt.Sprintf("Alignment length is not 184 (%d)", align.Length()))
	}
	if align.NbSequences() != 6 {
		t.Error(fmt.Sprintf("There are not 6 sequence in the alignment (%d)", align.NbSequences()))
	}

	align, err = NewParser(strings.NewReader(clustalstring2)).Parse()

	if err != nil {
		t.Error(err)
	}
	if align.Length() != 30 {
		t.Error(fmt.Sprintf("Alignment length is not 184 (%d)", align.Length()))
	}
	if align.NbSequences() != 6 {
		t.Error(fmt.Sprintf("There are not 6 sequence in the alignment (%d)", align.NbSequences()))
	}
}

func TestParseError1(t *testing.T) {
	_, err := NewParser(strings.NewReader(clustalstring3)).Parse()

	if err == nil {
		t.Error("There should be an error while reading this alignment")
	}
}

func TestParseError2(t *testing.T) {
	_, err := NewParser(strings.NewReader(clustalstring4)).Parse()

	if err == nil {
		t.Error("There should be an error while reading this alignment")
	}
}
