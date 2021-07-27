package draw

import (
	"bufio"
	"fmt"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io/fasta"
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
func (layout *bioJSLayout) DrawAlign(a align.Alignment) (err error) {
	var biojs string

	if biojs, err = biojsDepString(); err != nil {
		return
	}

	layout.writer.WriteString(`
<html>
<head>
  <script>
`)
	layout.writer.WriteString(biojs)
	layout.writer.WriteString(`
  </script>
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

	layout.writer.WriteString(`
  </pre>
  <div id="align" style="width:100%;height:100%;">
  </div>
    <script>
    // parsed array of the sequences
    var fastaseq = document.getElementById("fasta").innerText;
    var seqs =  msa.io.fasta.parse(fastaseq);
    var m = msa({
         el: align,
         seqs: seqs,
         vis: {
           conserv: true,
           overviewbox: false
         },
         // smaller menu for JSBin
         menu: "small",
         bootstrapMenu: true
    });
    m.g.zoomer.set("alignmentHeight", 500); // modifies the default height
    m.g.zoomer.set("labelNameLength", `)
	layout.writer.WriteString(fmt.Sprintf("%d);", a.MaxNameLength()*10))
	layout.writer.WriteString(`
    m.render();
   </script>
</body>
</html>
`)
	return err
}
