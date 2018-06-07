package draw

import (
	"bufio"

	"github.com/fredericlemoine/goalign/align"
	"github.com/fredericlemoine/goalign/io/fasta"
)

type bioJSLayout struct {
	writer *bufio.Writer
}

func NewBioJSLayout(writer *bufio.Writer) AlignLayout {
	return &bioJSLayout{writer}
}

/*
Draw the tree on the specific drawer. Does not close the file. The caller must do it.
*/
func (layout *bioJSLayout) DrawAlign(a align.Alignment) error {
	var err error = nil
	_, err = layout.writer.WriteString(`
<html>
<head>
  <script src="https://s3.eu-central-1.amazonaws.com/cdn.bio.sh/msa/latest/msa.min.gz.js"></script>
  <style media="screen" type="text/css">
    #cy {
    width: 100%;
    height: 100%;
    display: block;
    }
  </style>
</head>
<body>
  <pre id="fasta" style="display: none;">
`)
	layout.writer.WriteString(fasta.WriteAlignment(a))
	_, err = layout.writer.WriteString(`
  </pre>
  <div id="align" style="width:100%;height:100%;">
  </div>
    <script>
    // parsed array of the sequences
    var fastaseq = document.getElementById("fasta").innerText;
    var seqs =  msa.io.fasta.parse(fastaseq);
    var m = msa({
         el: align,
         seqs: seqs
    });
    m.g.zoomer.set("alignmentHeight", 500); // modifies the default height
    m.render();
   </script>
</body>
</html>
`)
	return err
}