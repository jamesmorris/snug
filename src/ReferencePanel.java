import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.swing.JPanel;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class ReferencePanel extends JPanel {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	ArrayList<SAMRecord> reads;
	int pixPerBase;
	Variant variant;
	int trackHeight;
	IndexedFastaSequenceFile reference;
	int start;
	int end;
	
	public ReferencePanel(int pixPerBase, Variant v, IndexedFastaSequenceFile ref, int start, int end) {
		super();
		this.pixPerBase = pixPerBase;
		this.variant = v;
		this.start = start;
		this.end = end;
		this.reference = ref;
	}

	public void paintComponent(Graphics g) {
		super.paintComponent(g);
		Graphics2D g2 = (Graphics2D) g;
		drawReference(g2);
	}

	private void drawReference(Graphics2D g2) {
		
		FontMetrics fm = g2.getFontMetrics();
		
		ReferenceSequence refSeq = reference.getSubsequenceAt(variant.getChr(), start, end);
		byte[] bases = refSeq.getBases();
		
		int ypos = 15;
		
		for(int i=0; i<bases.length; i++) {
			String base = String.valueOf((char) bases[i]);
			baseColour(g2, base);
			g2.drawString(base, i * pixPerBase, ypos);	
		}
		
		// draw a ruler
		int tickLength = 5;
		int basesPerTick = 20;
		ypos += 20;

		g2.setColor(Color.black);
		g2.setStroke(new BasicStroke(2));
		g2.drawLine(0, ypos, this.getWidth(), ypos);		
		g2.drawLine(1, ypos - 10, 1, ypos + 10);
		g2.drawLine(this.getWidth()-1, ypos - 10, this.getWidth()-1, ypos + 10);
		
		for(int i=0; i<bases.length; i++) {			
			if ((start + i) % basesPerTick == 0) {
				
				// if the genome position is divisable by by ten
				String bpString = Integer.toString(start + i);
				
				int wordWidth = fm.stringWidth(bpString);
				
				g2.drawLine(i * pixPerBase, ypos, i * pixPerBase, ypos + tickLength);
				g2.drawString(Integer.toString(start + i) + "bp", (i * pixPerBase) - (wordWidth / 2), ypos + 20);
			}			
		}
		
	}

	private void baseColour(Graphics2D g2, String base) {
		if (base.equals("A")) {
			g2.setColor(Color.green);
		} else if (base.equals("T")) {
			g2.setColor(Color.red);
		} else if (base.equals("C")) {
			g2.setColor(Color.blue);
		} else if (base.equals("G")) {
			g2.setColor(Color.orange);
		} else {
			g2.setColor(Color.black);
		}
	}

	public void setReference(IndexedFastaSequenceFile indexedFastaSequenceFile) {
		this.reference = indexedFastaSequenceFile;
	}

}
