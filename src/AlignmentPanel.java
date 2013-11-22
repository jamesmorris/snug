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
	private boolean mappingQuality;

	private Color vlightgrey = new Color(160, 160, 160);
	private Color lightgrey = new Color(192, 192, 192);
	private Color medgrey = new Color(224, 224, 224);
	private Color grey = Color.gray;
	private Color vlightgreen = new Color(204, 255, 204);
	private Color lightgreen = new Color(102, 255, 102);
	private Color medgreen = new Color(51, 255, 51);
	private Color green = Color.green;
	private Color vlightred = new Color(255, 204, 204);
	private Color lightred = new Color(255, 102, 102);
	private Color medred = new Color(255, 51, 51);
	private Color red = Color.red;
	private Color vlightorange = new Color(255, 229, 204);
	private Color lightorange = new Color(255, 204, 153);
	private Color medorange = new Color(255, 178, 102);
	private Color orange = new Color(255, 128, 0);
	private Color vlightblue = new Color(204, 229, 255);
	private Color lightblue = new Color(153, 204, 255);
	private Color medblue = new Color(102, 178, 255);
	private Color blue = Color.blue;

	ArrayList<BasePosition> baseCounts;

	public AlignmentPanel(int pixPerBase, int offset, ArrayList<SAMRecord> reads, Variant v, IndexedFastaSequenceFile ref, boolean collapse, boolean hide, boolean mappingQuality) {
		super();
		this.pixPerBase = pixPerBase;
		this.offset = offset;
		this.reads = reads;
		this.variant = v;
		this.reference = ref;
		this.collapse = collapse;
		this.hide = hide;
		this.mappingQuality = mappingQuality;

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
			drawAlignment(g2, 5);
		} else if (hide) {
			drawAlignment(g2, 15);
		} else {
			drawAlignment(g2, 15);
		}

	}

	private void drawAlignment(Graphics2D g2, int lineSpace) {
		if (reads != null) {
			int ypos = lineSpace;
			baseCounts = new ArrayList<BasePosition>();
			g2.setColor(Color.black);
			for (SAMRecord read : reads) {
				drawSequence(g2, read, ypos, offset);
				ypos += lineSpace;
			}
			// draw lines before and after to show the position of the variant
			g2.drawLine(offset - 2, 0, offset - 2, this.getHeight());
			g2.drawLine(offset + (pixPerBase - 2), 0, offset + (pixPerBase - 2), this.getHeight());
		}
	}

	private void drawSequence(final Graphics2D g2, final SAMRecord samRecord, int ypos, int offset) {

		int xpos = 0;
		int basePosition;
		final String readSeq = samRecord.getReadString();

		FontMetrics fm = g2.getFontMetrics();

		byte[] baseQualityScores = samRecord.getBaseQualities();
		int mappingQualityScore = samRecord.getMappingQuality();

		List<AlignmentBlock> blocks = samRecord.getAlignmentBlocks();

		// print the strand direction
		boolean strand = samRecord.getReadNegativeStrandFlag();

		if (strand) {
			samRecord.getAlignmentStart();
			int distance2Variant = samRecord.getAlignmentStart() - variant.getPos();
			basePosition = distance2Variant * pixPerBase;
			xpos = (offset + basePosition) - ((int) fm.getStringBounds("< < <", g2).getWidth() + 5);
			g2.setColor(Color.gray);
			if (collapse) {
				//
			} else if (hide) {
				//
			} else {
				g2.drawString("< < <", xpos, ypos);
			}
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

					if (mappingQuality) {
						refMappingQualityColour(g2, mappingQualityScore);
					} else {
						refBaseQualityColour(g2, baseQualityScores[readPos]);
					}

					if (collapse || hide) {
						g2.drawLine(xpos, ypos, xpos + pixPerBase, ypos);
					} else {
						g2.drawString(readBase, xpos, ypos);
					}
				} else {

					if (mappingQuality) {
						altMappingQualityColour(g2, readBase, mappingQualityScore);
					} else {
						altBaseQualityColour(g2, readBase, baseQualityScores[readPos]);
					}
					if (collapse) {
						g2.drawLine(xpos, ypos, xpos + pixPerBase, ypos);
					} else {
						g2.drawString(readBase, xpos, ypos);
					}
				}

				// update the base counts for the current position
				int position = xpos / pixPerBase;
				try {
					// try to increment the base count					
					if (strand) {
						baseCounts.get(position).addBase("-" + readBase);
					} else {
						baseCounts.get(position).addBase("+" + readBase);
					}
				} catch (IndexOutOfBoundsException e) {
					baseCounts.add(position, new BasePosition());
					if (strand) {
						baseCounts.get(position).addBase("-" + readBase);
					} else {
						baseCounts.get(position).addBase("+" + readBase);
					}
				}
			}
		}

		// print the strand direction
		if (!strand) {
			g2.setColor(Color.gray);
			if (collapse) {
				//
			} else if (hide) {
				//
			} else {
				g2.drawString("> > >", xpos + pixPerBase, ypos);
			}
		}

	}

	private void refMappingQualityColour(Graphics2D g2, int mappingQualityScore) {
		if (mappingQualityScore < 10)
			g2.setColor(vlightgrey);
		else if (mappingQualityScore < 20)
			g2.setColor(lightgrey);
		else if (mappingQualityScore < 30)
			g2.setColor(medgrey);
		else
			g2.setColor(grey);
	}

	private void refBaseQualityColour(Graphics2D g2, byte baseQuality) {
		if (baseQuality < 10)
			g2.setColor(vlightgrey);
		else if (baseQuality < 20)
			g2.setColor(lightgrey);
		else if (baseQuality < 30)
			g2.setColor(medgrey);
		else
			g2.setColor(grey);
	}

	private void altMappingQualityColour(Graphics2D g2, String base, int mappingQualityScore) {
		
		if (mappingQualityScore < 10){
			if (base.equals("A")) {			
				g2.setColor(vlightgreen);
			} else if (base.equals("T")) {
				g2.setColor(vlightred);
			} else if (base.equals("C")) {
				g2.setColor(vlightblue);
			} else if (base.equals("G")) {
				g2.setColor(vlightorange);
			} else {
				g2.setColor(Color.black);
			}
		} else if (mappingQualityScore < 20) {
			if (base.equals("A")) {			
				g2.setColor(lightgreen);
			} else if (base.equals("T")) {
				g2.setColor(lightred);
			} else if (base.equals("C")) {
				g2.setColor(lightblue);
			} else if (base.equals("G")) {
				g2.setColor(lightorange);
			} else {
				g2.setColor(Color.black);
			}
		} else if (mappingQualityScore < 30) {
			if (base.equals("A")) {			
				g2.setColor(medgreen);
			} else if (base.equals("T")) {
				g2.setColor(medred);
			} else if (base.equals("C")) {
				g2.setColor(medblue);
			} else if (base.equals("G")) {
				g2.setColor(medorange);
			} else {
				g2.setColor(Color.black);
			}
		} else {
			if (base.equals("A")) {			
				g2.setColor(green);
			} else if (base.equals("T")) {
				g2.setColor(red);
			} else if (base.equals("C")) {
				g2.setColor(blue);
			} else if (base.equals("G")) {
				g2.setColor(orange);
			} else {
				g2.setColor(Color.black);
			}
		}
	}

	private void altBaseQualityColour(Graphics2D g2, String base, byte baseQuality) {
		if (baseQuality < 10){
			if (base.equals("A")) {			
				g2.setColor(vlightgreen);
			} else if (base.equals("T")) {
				g2.setColor(vlightred);
			} else if (base.equals("C")) {
				g2.setColor(vlightblue);
			} else if (base.equals("G")) {
				g2.setColor(vlightorange);
			} else {
				g2.setColor(Color.black);
			}
		} else if (baseQuality < 20) {
			if (base.equals("A")) {			
				g2.setColor(lightgreen);
			} else if (base.equals("T")) {
				g2.setColor(lightred);
			} else if (base.equals("C")) {
				g2.setColor(lightblue);
			} else if (base.equals("G")) {
				g2.setColor(lightorange);
			} else {
				g2.setColor(Color.black);
			}
		} else if (baseQuality < 30) {
			if (base.equals("A")) {			
				g2.setColor(medgreen);
			} else if (base.equals("T")) {
				g2.setColor(medred);
			} else if (base.equals("C")) {
				g2.setColor(medblue);
			} else if (base.equals("G")) {
				g2.setColor(medorange);
			} else {
				g2.setColor(Color.black);
			}
		} else {
			if (base.equals("A")) {			
				g2.setColor(green);
			} else if (base.equals("T")) {
				g2.setColor(red);
			} else if (base.equals("C")) {
				g2.setColor(blue);
			} else if (base.equals("G")) {
				g2.setColor(orange);
			} else {
				g2.setColor(Color.black);
			}
		}
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
		return "<html><table>" +
		"<tr><td></td><td>+</td><td>-</td></tr>"+
		"<tr><td>A</td><td>" + bp.getAplus() + "</td><td>" + bp.getAminus() + "</td></tr>"+
		"<tr><td>T</td><td>" + bp.getTplus() + "</td><td>" + bp.getTminus() + "</td></tr>"+
		"<tr><td>G</td><td>" + bp.getGplus() + "</td><td>" + bp.getGminus() + "</td></tr>"+
		"<tr><td>C</td><td>" + bp.getCplus() + "</td><td>" + bp.getCminus() + "</td></tr>"+
		"</table></html>";
	}

}
