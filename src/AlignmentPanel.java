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
	IndexedFastaSequenceFile reference;
	// interface options
	private boolean collapse;
	private boolean hide;
	
	ArrayList<BasePosition> baseCounts;

	public AlignmentPanel(int pixPerBase, int offset, ArrayList<SAMRecord> reads, Variant v, IndexedFastaSequenceFile ref, boolean collapse, boolean hide) {
		super();
		this.pixPerBase = pixPerBase;
		this.offset = offset;
		this.reads = reads;
		this.variant = v;
		this.reference = ref;
		this.collapse = collapse;
		this.hide = hide;

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
		if (collapse) {
			drawCollapsedAlignment(g2);
		} else if (hide) {
			drawHiddenAlignment(g2);
		} else {
			drawAlignment(g2);
		}
	}

	private void drawHiddenAlignment(Graphics2D g2) {
		if (reads != null) {
			int ypos = 15;

			baseCounts = new ArrayList<BasePosition>();

			g2.setColor(Color.black);
			for (SAMRecord read : reads) {
				drawHiddenSequence(g2, read, ypos, offset);
				ypos += 15;
			}
			// draw lines before and after to show the position of the variant
			g2.drawLine(offset - 2, 0, offset - 2, this.getHeight());
			g2.drawLine(offset + (pixPerBase - 2), 0, offset + (pixPerBase - 2), this.getHeight());
		}
		
	}

	private void drawHiddenSequence(Graphics2D g2, SAMRecord samRecord, int ypos, int offset2) {
		int xpos = 0;
		int basePosition;
		final String readSeq = samRecord.getReadString();

		byte[] baseQualityScores = samRecord.getBaseQualities();
		
		// set the line width
		g2.setStroke(new BasicStroke(2));

		List<AlignmentBlock> blocks = samRecord.getAlignmentBlocks();

		for (int i = 0; i < blocks.size(); i++) {

			AlignmentBlock block = blocks.get(i);
			int blockStart = block.getReadStart();
			int distance2Variant = block.getReferenceStart() - variant.getPos();
			for (int j = 0; j < block.getLength(); j++) {

				int readPos = blockStart - 1 + j;
				int currentPosition = ((distance2Variant + j) + variant.getPos());
				ReferenceSequence refSeq = reference.getSubsequenceAt(variant.getChr(), currentPosition, currentPosition);
				byte[] bases = refSeq.getBases();
				String refBase = String.valueOf((char) bases[0]);
				String readBase = readSeq.substring(readPos, readPos + 1);
				basePosition = (distance2Variant + j) * pixPerBase;
				xpos = offset + basePosition;

				// check read base against reference base - colour if different
				if (readBase.matches(refBase)) {
					refBaseQualityColour(g2, baseQualityScores[readPos]);
					g2.drawLine(xpos, ypos, xpos + pixPerBase, ypos);
				} else {
					altBaseQualityColour(g2, readBase, baseQualityScores[readPos]);
					g2.drawString(readBase, xpos, ypos);
				}

				// update the base counts for the current position
				int position = xpos / pixPerBase;

				try {
					// try to increment the base count
					baseCounts.get(position).addBase(readBase);
				} catch (IndexOutOfBoundsException e) {
					baseCounts.add(position, new BasePosition());
					baseCounts.get(position).addBase(readBase);
				}
			}
		}
		
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
			g2.drawLine(offset - 2, 0, offset - 2, this.getHeight());
			g2.drawLine(offset + (pixPerBase - 2), 0, offset + (pixPerBase - 2), this.getHeight());
		}
	}

	private void drawCollapsedAlignment(Graphics2D g2) {
		if (reads != null) {
			int ypos = 5;

			baseCounts = new ArrayList<BasePosition>();

			g2.setColor(Color.black);
			for (SAMRecord read : reads) {
				drawCollapsedSequence(g2, read, ypos, offset);
				ypos += 5;
			}
			// draw lines before and after to show the position of the variant
			g2.drawLine(offset - 2, 0, offset - 2, this.getHeight());
			g2.drawLine(offset + (pixPerBase - 2), 0, offset + (pixPerBase - 2), this.getHeight());
		}
	}

	private void drawCollapsedSequence(Graphics2D g2, SAMRecord samRecord, int ypos, int offset2) {

		int xpos = 0;
		int basePosition;
		final String readSeq = samRecord.getReadString();

		byte[] baseQualityScores = samRecord.getBaseQualities();
		
		// set the line width
		g2.setStroke(new BasicStroke(2));

		List<AlignmentBlock> blocks = samRecord.getAlignmentBlocks();

		for (int i = 0; i < blocks.size(); i++) {

			AlignmentBlock block = blocks.get(i);
			int blockStart = block.getReadStart();
			int distance2Variant = block.getReferenceStart() - variant.getPos();
			for (int j = 0; j < block.getLength(); j++) {

				int readPos = blockStart - 1 + j;
				int currentPosition = ((distance2Variant + j) + variant.getPos());
				ReferenceSequence refSeq = reference.getSubsequenceAt(variant.getChr(), currentPosition, currentPosition);
				byte[] bases = refSeq.getBases();
				String refBase = String.valueOf((char) bases[0]);
				String readBase = readSeq.substring(readPos, readPos + 1);
				basePosition = (distance2Variant + j) * pixPerBase;
				xpos = offset + basePosition;

				// check read base against reference base - colour if different
				if (readBase.matches(refBase)) {
					refBaseQualityColour(g2, baseQualityScores[readPos]);
				} else {
					altBaseQualityColour(g2, readBase, baseQualityScores[readPos]);
				}

				// update the base counts for the current position
				int position = xpos / pixPerBase;

				try {
					// try to increment the base count
					baseCounts.get(position).addBase(readBase);
				} catch (IndexOutOfBoundsException e) {
					baseCounts.add(position, new BasePosition());
					baseCounts.get(position).addBase(readBase);
				}
				g2.drawLine(xpos, ypos, xpos + pixPerBase, ypos);
			}
		}

	}

	private void drawSequence(final Graphics2D g2, final SAMRecord samRecord, int ypos, int offset) {

		int xpos = 0;
		int basePosition;
		final String readSeq = samRecord.getReadString();

		FontMetrics fm = g2.getFontMetrics();

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

		for (int i = 0; i < blocks.size(); i++) {

			AlignmentBlock block = blocks.get(i);
			int blockStart = block.getReadStart();
			int distance2Variant = block.getReferenceStart() - variant.getPos();
			for (int j = 0; j < block.getLength(); j++) {

				int readPos = blockStart - 1 + j;
				int currentPosition = ((distance2Variant + j) + variant.getPos());
				ReferenceSequence refSeq = reference.getSubsequenceAt(variant.getChr(), currentPosition, currentPosition);
				byte[] bases = refSeq.getBases();
				String refBase = String.valueOf((char) bases[0]);
				String readBase = readSeq.substring(readPos, readPos + 1);
				basePosition = (distance2Variant + j) * pixPerBase;
				xpos = offset + basePosition;

				// check read base against reference base - colour if different
				if (readBase.matches(refBase)) {
					refBaseQualityColour(g2, baseQualityScores[readPos]);
				} else {
					altBaseQualityColour(g2, readBase, baseQualityScores[readPos]);
				}
				
				// update the base counts for the current position
				int position = xpos / pixPerBase;
				try {
					// try to increment the base count
					baseCounts.get(position).addBase(readBase);
				} catch (IndexOutOfBoundsException e) {
					baseCounts.add(position, new BasePosition());
					baseCounts.get(position).addBase(readBase);
				}
				g2.drawString(readBase, xpos, ypos);

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
	
	private void refBaseQualityColour(Graphics2D g2, byte baseQuality) {
		Color vlightgrey = new Color(160, 160, 160);
		Color lightgrey = new Color(192, 192, 192);
		Color medgrey = new Color(224, 224, 224);
		Color grey = Color.gray;
		
		if (baseQuality < 10)
			g2.setColor(vlightgrey);
		else if (baseQuality < 20)
			g2.setColor(lightgrey);
		else if (baseQuality < 30)
			g2.setColor(medgrey);
		else
			g2.setColor(grey);
	}
	
	private void altBaseQualityColour(Graphics2D g2, String base, byte baseQuality) {
		if (base.equals("A")) {			
			setAColourByBaseQuality(g2, baseQuality);
		} else if (base.equals("T")) {
			setTColourByBaseQuality(g2, baseQuality);
		} else if (base.equals("C")) {
			setCColourByBaseQuality(g2, baseQuality);
		} else if (base.equals("G")) {
			setGColourByBaseQuality(g2, baseQuality);
		} else {
			g2.setColor(Color.black);
		}
	}
	
	private void setAColourByBaseQuality(Graphics2D g2, byte baseQuality) {
		Color vlightgreen = new Color(204, 255, 204);
		Color lightgreen = new Color(102, 255, 102);
		Color medgreen = new Color(51, 255, 51);
		Color green = Color.green;
		
		if (baseQuality < 10)
			g2.setColor(vlightgreen);
		else if (baseQuality < 20)
			g2.setColor(lightgreen);
		else if (baseQuality < 30)
			g2.setColor(medgreen);
		else
			g2.setColor(green);
	}

	private void setTColourByBaseQuality(Graphics2D g2, byte baseQuality) {
		Color vlightred = new Color(255, 204, 204);
		Color lightred = new Color(255, 102, 102);
		Color medred = new Color(255, 51, 51);
		Color red = Color.red;
		
		if (baseQuality < 10)
			g2.setColor(vlightred);
		else if (baseQuality < 20)
			g2.setColor(lightred);
		else if (baseQuality < 30)
			g2.setColor(medred);
		else
			g2.setColor(red);
	}

	private void setGColourByBaseQuality(Graphics2D g2, byte baseQuality) {
		Color vlightorange = new Color(255, 229, 204);
		Color lightorange = new Color(255, 204, 153);
		Color medorange = new Color(255, 178, 102);
		Color orange = new Color(255, 128, 0);
		
		if (baseQuality < 10)
			g2.setColor(vlightorange);
		else if (baseQuality < 20)
			g2.setColor(lightorange);
		else if (baseQuality < 30)
			g2.setColor(medorange);
		else
			g2.setColor(orange);
	}

	private void setCColourByBaseQuality(Graphics2D g2, byte baseQuality) {
		Color vlightblue = new Color(204, 229, 255);
		Color lightblue = new Color(153, 204, 255);
		Color medblue = new Color(102, 178, 255);
		Color blue = Color.blue;
		
		if (baseQuality < 10)
			g2.setColor(vlightblue);
		else if (baseQuality < 20)
			g2.setColor(lightblue);
		else if (baseQuality < 30)
			g2.setColor(medblue);
		else
			g2.setColor(blue);
	}

	public void setReads(ArrayList<SAMRecord> reads) {
		this.reads = reads;
	}

	public void setReference(IndexedFastaSequenceFile indexedFastaSequenceFile) {
		this.reference = indexedFastaSequenceFile;
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
