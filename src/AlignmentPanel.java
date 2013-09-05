import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.event.MouseEvent;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.swing.JPanel;
import javax.swing.ToolTipManager;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class AlignmentPanel extends JPanel {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	ArrayList<SAMRecord> reads;
	int pixPerBase;
	int offset;
	Variant variant;
	int trackHeight;
	IndexedFastaSequenceFile reference;
	// interface options
	private boolean hideBases;
	
	ArrayList<BasePosition> baseCounts;
	
	public AlignmentPanel(int pixPerBase, int offset, ArrayList<SAMRecord> reads, Variant v, int height, IndexedFastaSequenceFile ref) {
		super();
		this.pixPerBase = pixPerBase;
		this.offset = offset;
		this.reads = reads;
		this.variant = v;
		this.trackHeight = height;
		this.reference = ref;
		
		ToolTipManager.sharedInstance().registerComponent(this);
		ToolTipManager.sharedInstance().setDismissDelay(10000);
		ToolTipManager.sharedInstance().setInitialDelay(0);
		ToolTipManager.sharedInstance().setReshowDelay(0);
			
	}

	public AlignmentPanel() {
		// TODO Auto-generated constructor stub
	}

	public void paintComponent(Graphics g) {
		super.paintComponent(g);
		Graphics2D g2 = (Graphics2D) g;
		drawAlignment(g2);
	}
	
	private void drawAlignment(Graphics2D g2) {
    	if (reads != null) {    		
    		int ypos = 15;
    		
    		baseCounts = new ArrayList<BasePosition>();
    		
    		g2.setColor(Color.black);
            for (SAMRecord read : reads) {            	
            	drawSequence(g2, read, ypos, offset);
            	ypos += 15;
    		}
    		// draw lines before and after to show the position of the variant
            g2.drawLine(offset-2, 0, offset-2, trackHeight);
    		g2.drawLine(offset + (pixPerBase-2), 0, offset + (pixPerBase-2), trackHeight);
    	}
    }
	
	private void drawSequence(final Graphics2D g2, final SAMRecord samRecord, int ypos, int offset) {
		
		int xpos = 0;
		int basePosition;
		final String readSeq = samRecord.getReadString();
		
		FontMetrics fm = g2.getFontMetrics();
		Rectangle2D rect = fm.getStringBounds("A", g2);

		byte[] baseQualityScores = samRecord.getBaseQualities();

		List<AlignmentBlock> blocks = samRecord.getAlignmentBlocks();
		
		// print the strand direction		
		boolean strand = samRecord.getReadNegativeStrandFlag();
		
		if (strand) {
			samRecord.getAlignmentStart();
			int distance2Variant = samRecord.getAlignmentStart() - variant.getPos();			
			basePosition = distance2Variant * pixPerBase;			
			xpos = (offset + basePosition) - ((int) fm.getStringBounds("< < <", g2).getWidth() + 5);
			g2.setColor(Color.gray);
			g2.drawString("< < <", xpos, ypos);			
		}
		
		for(int i=0; i<blocks.size(); i++) {
			
			AlignmentBlock block = blocks.get(i);
			int blockStart = block.getReadStart();
			int distance2Variant = block.getReferenceStart() - variant.getPos();			
			for(int j=0; j<block.getLength(); j++) {
			
				int readPos = blockStart-1+j;				
				int currentPosition = ((distance2Variant + j) + variant.getPos());				
				ReferenceSequence refSeq = reference.getSubsequenceAt(variant.getChr(), currentPosition, currentPosition);
				byte[] bases = refSeq.getBases();				
				String refBase = String.valueOf((char) bases[0]);
				String readBase = readSeq.substring(readPos, readPos+1);
				basePosition = (distance2Variant + j) * pixPerBase;				
				xpos = offset + basePosition;
				
				// shade the background of the base on quality
				setColourByBaseQuality(g2, baseQualityScores[readPos]);
				rect = fm.getStringBounds(readBase, g2);
				g2.fillRect(xpos, ypos - fm.getAscent(), (int) rect.getWidth(), (int) rect.getHeight());
				
				// check read base against reference base - colour if different
				if (readBase.matches(refBase)) {
					g2.setColor(Color.gray);					
				} else {
					baseColour(g2, readBase);
				}				
				
				// update the base counts for the current position
				int position = xpos / pixPerBase;
				
				// first test if the base is in view
				if (xpos >= 0) {
					try {
						baseCounts.get(position).addBase(readBase);		
					} catch (IndexOutOfBoundsException e) {
						baseCounts.add(position, new BasePosition());
						baseCounts.get(position).addBase(readBase);
					}
					
					g2.drawString(readBase, xpos, ypos);	
				}
			}			
		}
	
		// print the strand direction
		if (!strand) {
			g2.setColor(Color.gray);
			g2.drawString("> > >", xpos + pixPerBase, ypos);
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

	private void setColourByBaseQuality(Graphics2D g2, byte baseQuality) {
	    if (baseQuality < 10)
	      g2.setColor(new Color(160, 160, 160));
	    else if (baseQuality < 20)
	      g2.setColor(new Color(192, 192, 192));
	    else if (baseQuality < 30)
	      g2.setColor(new Color(224, 224, 224));
	    else
	    	g2.setColor(Color.white);
	}
	
	public void setReads(ArrayList<SAMRecord> reads) {
		this.reads = reads;
	}

	public void setReference(IndexedFastaSequenceFile indexedFastaSequenceFile) {
		this.reference = indexedFastaSequenceFile;
	}

	public void setHideBases(boolean b) {
		this.hideBases = b;		
	}

    @Override
    public String getToolTipText(MouseEvent e) {
		int basePosition = e.getX() / pixPerBase;
		return getBaseCounts(basePosition);
    }

	private String getBaseCounts(int basePosition) {
		BasePosition bp = baseCounts.get(basePosition);
		return "A:" + bp.getA() + " T" + bp.getT() + " C:" + bp.getC() + " G:" + bp.getG();
	}
	
}
