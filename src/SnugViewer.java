import javax.imageio.ImageIO;
import javax.swing.*;
import javax.swing.border.Border;
import javax.swing.border.EmptyBorder;
import javax.swing.border.LineBorder;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;

import java.awt.event.*;
import java.awt.geom.Ellipse2D;
import java.awt.image.BufferedImage;
import java.awt.*;
import java.io.*;
import java.text.NumberFormat;
import java.util.Arrays;
import java.util.Collections;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;
import java.util.Locale;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.concurrent.TimeUnit;


public class SnugViewer extends JFrame implements ActionListener, ComponentListener {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1122484803356474269L;
	String plottedSNP = null;
	private File output;
	private ArrayList<Variant> variantList;
	private Variant currentVariantinList;
	private Variant currentVariant = null;
	private int currentVariantIndex;
	private Integer displayVariantIndex = null;
	private boolean backVariant = false;
	private boolean bamFile = false;
	private boolean variantFile = false;
	private boolean endOfList;
	private HashMap<Variant, Integer> listScores = new HashMap<Variant, Integer>();
	JFileChooser jfc;
	private JPanel scorePanel;
	private JPanel messagePanel;
	private JPanel controlsPanel;
	private JLabel message;
	private JButton yesButton;
	private JButton maybeButton;
	private JButton noButton;
	private JButton backButton;
	private JMenu fileMenu;
	private JMenuItem loadFiles;
	private JMenuItem showLogItem;
	public static LoggingDialog ld;
	DataConnectionDialog dcd;
	
	JMenuBar mb;
	
	// untested below
	private JPanel trackContainer;
	
	private JCheckBoxMenuItem hideBases;
	
	JPanel contentPanel;

	SAMFileReader bamReader;
	private ArrayList<SAMFileReader> bams;

	private boolean backSNP;
	
	private IndexedFastaSequenceFile referenceFile;
	

	int pixPerBase = 15;
	

	public static void main(String[] args) {

		new SnugViewer();

	}

	SnugViewer() {
		super("Snug...");

		dcd = new DataConnectionDialog(this);
		mb = new JMenuBar();

		int menumask = Toolkit.getDefaultToolkit().getMenuShortcutKeyMask();

		fileMenu = new JMenu("File");
		loadFiles = new JMenuItem("Load Files");
		loadFiles.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_B, menumask));
		loadFiles.addActionListener(this);
		fileMenu.add(loadFiles);
		
		hideBases = new JCheckBoxMenuItem("Hide Bases");
		hideBases.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_F, menumask));
		hideBases.addActionListener(this);
		hideBases.setEnabled(false);
		fileMenu.add(hideBases);		
		
		mb.add(fileMenu);

		if (!(System.getProperty("os.name").toLowerCase().contains("mac"))) {
			JMenuItem quitItem = new JMenuItem("Quit");
			quitItem.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_Q,menumask));
			quitItem.addActionListener(this);
			fileMenu.add(quitItem);
		}

		JMenu logMenu = new JMenu("Log");
		showLogItem = new JMenuItem("Show log");
		showLogItem.addActionListener(this);
		logMenu.add(showLogItem);
		mb.add(logMenu);

		setJMenuBar(mb);

		controlsPanel = new JPanel();

		scorePanel = new JPanel();
		scorePanel.add(new JLabel("Approve?"));

		yesButton = new JButton("Yes");
		scorePanel.registerKeyboardAction(this, "Yes", KeyStroke.getKeyStroke(KeyEvent.VK_Y, 0), JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT);
		yesButton.addActionListener(this);
		yesButton.setEnabled(false);
		scorePanel.add(yesButton);

		maybeButton = new JButton("Maybe");
		scorePanel.registerKeyboardAction(this, "Maybe", KeyStroke.getKeyStroke(KeyEvent.VK_M, 0), JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT);
		maybeButton.addActionListener(this);
		maybeButton.setEnabled(false);
		scorePanel.add(maybeButton);

		noButton = new JButton("No");
		scorePanel.registerKeyboardAction(this, "No", KeyStroke.getKeyStroke(KeyEvent.VK_N, 0), JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT);
		noButton.addActionListener(this);
		noButton.setEnabled(false);
		scorePanel.add(noButton);

		backButton = new JButton("Back");
		scorePanel.registerKeyboardAction(this, "Back", KeyStroke.getKeyStroke(KeyEvent.VK_B, 0), JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT);
		backButton.addActionListener(this);
		backButton.setEnabled(false);
		scorePanel.add(backButton);

		messagePanel = new JPanel();
		message = new JLabel("");
		message.setEnabled(false);
		messagePanel.add(message);
		message.setVisible(false);

		JPanel rightPanel = new JPanel();
		rightPanel.add(scorePanel);
		rightPanel.add(messagePanel);
		rightPanel.setLayout(new BoxLayout(rightPanel, BoxLayout.Y_AXIS));
		controlsPanel.add(rightPanel);

		controlsPanel.setMaximumSize(new Dimension(2000, (int) controlsPanel.getPreferredSize().getHeight()));
		controlsPanel.setMinimumSize(new Dimension(10, (int) controlsPanel.getPreferredSize().getHeight()));

		trackContainer = new JPanel();
		trackContainer.setPreferredSize(new Dimension(900, 600));
		trackContainer.setBackground(Color.white);
		
		
		contentPanel = new JPanel();
		contentPanel.setLayout(new BoxLayout(contentPanel, BoxLayout.Y_AXIS));
		contentPanel.add(trackContainer);
		contentPanel.add(controlsPanel);
		contentPanel.addComponentListener(this);
		

		addWindowListener(new WindowAdapter() {

			public void windowClosing(WindowEvent e) {
				System.exit(0);
			}
		});

		this.setContentPane(contentPanel);
		this.pack();
		this.setVisible(true);

		ld = new LoggingDialog(this);
		ld.pack();
	}

	public void actionPerformed(ActionEvent actionEvent) {

		String noCommand = new String("No");
		String maybeCommand = new String("Maybe");
		String yesCommand = new String("Yes");
		String backCommand = new String("Back");
		String openFilesCommand = new String("Load Files");
		String showLogCommand = new String("Show log");
		String quitCommand = new String("Quit");
		String hideBasesCommand = "Hide Bases";

		String command = actionEvent.getActionCommand();
		if (command.equals(noCommand)) {
			noButton.requestFocusInWindow();
			try {
				recordVerdict(-1);
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		} else if (command.equals(maybeCommand)) {
			maybeButton.requestFocusInWindow();
			try {
				recordVerdict(0);
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		} else if (command.equals(yesCommand)) {
			yesButton.requestFocusInWindow();
			try {
				recordVerdict(1);
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		} else if (command.equals(backCommand)) {
			setBackSNP(true);
			displayVariantIndex--;
			if (displayVariantIndex == 0) {
				backButton.setEnabled(false);
			}
			int lastCall = listScores.get(variantList.get(displayVariantIndex));
			if (lastCall == 1) {
				yesButton.requestFocusInWindow();
			} else if (lastCall == 0) {
				maybeButton.requestFocusInWindow();
			} else {
				noButton.requestFocusInWindow();
			}
			drawReads();
		} else if (command.equals(openFilesCommand)) {

			bams = new ArrayList<SAMFileReader>();
			
			dcd.pack();
			dcd.setVisible(true);

			// check all the dialog fields are valid paths to data
			
			String bam = dcd.getBam();
			
			// if the bam path ends .bam just open the single file			
			if (bam.endsWith(".bam")) {
				System.out.println(bam);
				File bamFile = new File(bam);
				if (bamFile.exists() && checkBamIndex(bam)) {
					bams.add(new SAMFileReader(bamFile));
					bamReader = new SAMFileReader(bamFile);
				}
			
			// if the bam path ends .bams the file is assumed to contain a list of bam files				
			} else if (bam.endsWith(".bams")) {
				
				File bamFiles = new File(bam);
				String absolutePath = bamFiles.getAbsolutePath();
				String bamDir = absolutePath.substring(0,absolutePath.lastIndexOf(File.separator));
				
				// open the bams file and loop through all the bam files
				BufferedReader listReader = null;
				try {
					listReader = new BufferedReader(new FileReader(bam));
				} catch (FileNotFoundException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				
				String currentline;
				
				try {
					while ((currentline = listReader.readLine()) != null) {
						//check that the line is not blank
						if (currentline == null || currentline.trim().equals("")){
							continue;
						} else {
							String b = bamDir +File.separator + currentline.trim();
							File bamFile = new File(b);
							if (bamFile.exists() && checkBamIndex(b)) {
								bams.add(new SAMFileReader(bamFile));
							}
						}
					}
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}

			currentVariantIndex = 0;
			displayVariantIndex = currentVariantIndex;
			
			// load variants file
			String var = dcd.getVar();
			VariantListFile variantListFile = new VariantListFile(var);
			variantList = variantListFile.getVariants();
			ld.log(variantList.size() + " variants loaded");
			output = checkOverwriteFile(new File(var + ".scores"));
			
			
			// load the reference file			
			try {
				referenceFile = new IndexedFastaSequenceFile(new File(dcd.getRef()));
				ld.log("Reference file loaded");
			} catch (FileNotFoundException e) {
				ld.log("Reference file not found");
			}
			
			yesButton.setEnabled(true);
			noButton.setEnabled(true);
			maybeButton.setEnabled(true);
			yesButton.requestFocusInWindow();
			drawReads();
//			this.pack();
			
				
		} else if (command.equals(quitCommand)) {
			System.exit(0);
		} else if (command.equals(showLogCommand)) {
			ld.setVisible(true);
		} else if (command.equals(hideBasesCommand)) {
            // turn filtering on/off
            if (hideBases.isSelected()) {
                hideBases();
            } else {
                showBases();
            }
		}
		
	}

	private void showBases() {
		// TODO Auto-generated method stub
		
	}

	private void hideBases() {
		// TODO Auto-generated method stub
		
	}

	private Boolean checkBamIndex(String bam) {
		
		File bamIndexFile = new File(bam + ".bai");

		if (!bamIndexFile.exists()) {
			String ls = System.getProperty("line.separator");
			String msg = "BAM index file is missing. The BAM file needs to be sorted and indexed"
					+ ls
					+ "This can be done using samtools (http://samtools.sf.net/):"
					+ ls
					+ ls
					+ "samtools sort <in.bam> <out.prefix>"
					+ ls
					+ "samtools index <sorted.bam>";
			ld.log(msg);
			return Boolean.FALSE;
		} else {
			ld.log("BAM loaded: " + bam);
			return Boolean.TRUE;	
		}
		
	}

	private void drawReads() {
		Variant v = variantList.get(displayVariantIndex);
		this.setTitle(v.getName());
		
		trackContainer.removeAll();

		// calculate the height of each panel		
		int refPanelHeight = 60;
		int border = 40;
		int height = this.getHeight() - (controlsPanel.getHeight() + mb.getHeight()) - refPanelHeight - border;
		int trackHeight = height / bams.size();
		int width = this.getWidth() - border;
		int visableBases = width / pixPerBase;
		
		// assume the middle base is the variant - how many bases are either side of it
		int basesOneSide = (visableBases -1) / 2;
		
		// distance in pix to the variant site
		int offset = basesOneSide * pixPerBase;
				
		for (SAMFileReader sfr : bams) {
			
			ArrayList<SAMRecord> reads = new ArrayList<SAMRecord>();		
			SAMRecordIterator readsItr = sfr.queryOverlapping(v.getChr(), v.getPos(), v.getPos());
			
			while (readsItr.hasNext()) {	
	        	reads.add(readsItr.next());  	
	        }
	        readsItr.close();
	        
	        AlignmentPanel ap = new AlignmentPanel(pixPerBase, offset, reads, v, height, referenceFile);
	        
	        ap.setPreferredSize(new Dimension(width, trackHeight));
			setBackground(Color.white);
			ap.setBorder(BorderFactory.createLineBorder(Color.black));	        
			trackContainer.add(ap);

		}
		
		// draw the reference panel
		int refStart = v.getPos() - basesOneSide;
		int refEnd = v.getPos() + basesOneSide;

		ReferencePanel rp = new ReferencePanel(pixPerBase, v, referenceFile, refStart, refEnd);
		rp.setPreferredSize(new Dimension(width, refPanelHeight));
		trackContainer.add(rp);
		
		this.setVisible(true);
		
	}

	private void disableAllActions() {
		// disable all actions while loading data
		yesButton.setEnabled(false);
		noButton.setEnabled(false);
		maybeButton.setEnabled(false);
		backButton.setEnabled(false);
		scorePanel.setEnabled(false);
	}

	private void printScores() throws FileNotFoundException {
		ListIterator<Variant> var_i = variantList.listIterator();
		PrintStream ps = new PrintStream(new FileOutputStream(output));

		while (var_i.hasNext()) {
			Variant var = (Variant) var_i.next();
			if (listScores.get(var) != null) {
				ps.println(var.name + "\t" + listScores.get(var));
			}
		}
		ps.close();
	}

	private void activeScorePanel(boolean score) {
		if (score) {
			yesButton.setEnabled(true);
			noButton.setEnabled(true);
			maybeButton.setEnabled(true);
			scorePanel.setEnabled(true);
			if (currentVariantIndex == 0 || displayVariantIndex == 0) {
				backButton.setEnabled(false);
			} else {
				backButton.setEnabled(true);
			}
		} else {
			yesButton.setEnabled(false);
			noButton.setEnabled(false);
			maybeButton.setEnabled(false);
			backButton.setEnabled(false);
			scorePanel.setEnabled(false);
		}
	}

	private void setBackSNP(boolean back) {
		backSNP = back;
	}

	private boolean isBackSNP() {
		return backSNP;
	}

	private File checkOverwriteFile(File file) {

		if (file.exists()) {
			int n = JOptionPane
					.showConfirmDialog(
							this.getContentPane(),
							"The file "
									+ file.getName()
									+ " already exists\n would you like to overwrite this file?",
							"Overwrite file?", JOptionPane.YES_NO_OPTION,
							JOptionPane.QUESTION_MESSAGE);
			// n 0 = yes 1 = no
			if (n == 1) {
				if (jfc.showSaveDialog(this) == JFileChooser.APPROVE_OPTION) {
					file = new File(jfc.getSelectedFile().getAbsolutePath());
				} else {
					file = null;
				}
			}
		}
		return file;
	}

	private void recordVerdict(int v) throws FileNotFoundException {
		if (isBackSNP()) {
			listScores.put(variantList.get(displayVariantIndex), v);
			printScores();
			setBackSNP(false);
			if (displayVariantIndex == 0) {
				backButton.setEnabled(true);
			}
			displayVariantIndex = currentVariantIndex;
			drawReads();
		} else {
			listScores.put(variantList.get(currentVariantIndex), v);
			printScores();
			if (currentVariantIndex < (variantList.size() - 1)) {
				currentVariantIndex++;
				displayVariantIndex = currentVariantIndex;
				backButton.setEnabled(true);
				drawReads();
			} else {
				drawReads();
				activeScorePanel(false);
				int n = JOptionPane.showConfirmDialog(this.getContentPane(),
						"Would you like to finish scoring the current list?",
						"Finish scoring list?", JOptionPane.YES_NO_OPTION,
						JOptionPane.QUESTION_MESSAGE);
				// n 0 = yes 1 = no
				if (n == 0) {
					variantList.clear();
				} else {
					activeScorePanel(true);
					drawReads();
				}
			}
		}
	}

	private void printMessage(String string) {
		message.setText(string);
		message.setVisible(true);
	}

	public void refreshPlot() {
		if (displayVariantIndex != null) {
			drawReads();
		}
	}

	@Override
	public void componentHidden(ComponentEvent arg0) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void componentMoved(ComponentEvent arg0) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void componentResized(ComponentEvent arg0) {
		if (displayVariantIndex != null) {
			drawReads();
		}
	}

	@Override
	public void componentShown(ComponentEvent arg0) {
		// TODO Auto-generated method stub
		
	}
}
