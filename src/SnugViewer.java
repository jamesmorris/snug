import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.KeyStroke;

import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;
import java.awt.event.KeyEvent;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionListener;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;

import java.sql.Timestamp;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.HashMap;
import java.util.ListIterator;
import java.util.ArrayList;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.picard.reference.IndexedFastaSequenceFile;


public class SnugViewer extends JFrame implements ActionListener, ComponentListener {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1122484803356474269L;
	String plottedSNP = null;
	private File output;
	private ArrayList<Variant> variantList;
	private int currentVariantIndex;
	private Integer displayVariantIndex = null;
	private HashMap<Variant, Integer> listScores = new HashMap<Variant, Integer>();
	private JPanel scorePanel;
	private JScrollPane messagePanel;
	private JTextArea messageText;
	private JPanel controlsPanel;
	private JButton yesButton;
	private JButton maybeButton;
	private JButton noButton;
	private JButton backButton;
	private JMenu fileMenu;
	private JMenuItem loadFiles;
	DataConnectionDialog dcd;
	
	JMenuBar mb;
	
	private JPanel trackContainer;
	private JPanel nameContainer;
	private String refName;
	JPanel contentPanel;

	private ArrayList<SAMFileReader> bams;
	private ArrayList<String> bamNames;

	private boolean backSNP;
	
	private IndexedFastaSequenceFile referenceFile;
	
	// set the space given to each base as it is drawn
	int pixPerBase = 15;

	public static void main(String[] args) {
		new SnugViewer();
	}

	SnugViewer() {
		super("Snug...");

//		addMouseMotionListener(this);
		
		dcd = new DataConnectionDialog(this);
		mb = new JMenuBar();

		int menumask = Toolkit.getDefaultToolkit().getMenuShortcutKeyMask();

		fileMenu = new JMenu("File");
		loadFiles = new JMenuItem("Load Files");
		loadFiles.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_B, menumask));
		loadFiles.addActionListener(this);
		fileMenu.add(loadFiles);		
		
		mb.add(fileMenu);

		if (!(System.getProperty("os.name").toLowerCase().contains("mac"))) {
			JMenuItem quitItem = new JMenuItem("Quit");
			quitItem.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_Q,menumask));
			quitItem.addActionListener(this);
			fileMenu.add(quitItem);
		}

		setJMenuBar(mb);

		yesButton = new JButton("Yes");
		yesButton.addActionListener(this);
		yesButton.setEnabled(false);		

		maybeButton = new JButton("Maybe");
		maybeButton.addActionListener(this);
		maybeButton.setEnabled(false);
		
		noButton = new JButton("No");
		noButton.addActionListener(this);
		noButton.setEnabled(false);
		
		backButton = new JButton("Back");		
		backButton.addActionListener(this);
		backButton.setEnabled(false);
		
		scorePanel = new JPanel();
		scorePanel.add(new JLabel("Approve?"));
		scorePanel.add(yesButton);
		scorePanel.registerKeyboardAction(this, "Yes", KeyStroke.getKeyStroke(KeyEvent.VK_Y, 0), JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT);
		scorePanel.add(noButton);
		scorePanel.registerKeyboardAction(this, "No", KeyStroke.getKeyStroke(KeyEvent.VK_N, 0), JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT);		
		scorePanel.add(maybeButton);
		scorePanel.registerKeyboardAction(this, "Maybe", KeyStroke.getKeyStroke(KeyEvent.VK_M, 0), JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT);		
		scorePanel.add(backButton);
		scorePanel.registerKeyboardAction(this, "Back", KeyStroke.getKeyStroke(KeyEvent.VK_B, 0), JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT);
		
		messageText = new JTextArea(4, 40);
		messageText.setEditable(false);
		messagePanel = new JScrollPane(messageText);
		
		
		controlsPanel = new JPanel();
		controlsPanel.add(scorePanel);
		controlsPanel.add(messagePanel);
		controlsPanel.setLayout(new BoxLayout(controlsPanel, BoxLayout.X_AXIS));
		controlsPanel.setMaximumSize(new Dimension(900, 200));

		trackContainer = new JPanel();
		trackContainer.setPreferredSize(new Dimension(900, 600));
		trackContainer.setBackground(Color.white);
		trackContainer.setLayout(new BoxLayout(trackContainer, BoxLayout.Y_AXIS));
		
		nameContainer = new JPanel();
		nameContainer.setPreferredSize(new Dimension(60, 600));
		nameContainer.setBackground(Color.white);
		nameContainer.setLayout(new BoxLayout(nameContainer, BoxLayout.Y_AXIS));
		
		
		JPanel mainPanel = new JPanel();
		mainPanel.setLayout(new BoxLayout(mainPanel, BoxLayout.X_AXIS));
		mainPanel.add(nameContainer);
		mainPanel.add(trackContainer);
		
		contentPanel = new JPanel();
		contentPanel.setLayout(new BoxLayout(contentPanel, BoxLayout.Y_AXIS));
		contentPanel.add(mainPanel);
		contentPanel.add(controlsPanel);
		contentPanel.addComponentListener(this);

		controlsPanel.setVisible(false);
		
		this.setContentPane(contentPanel);
		this.pack();
		this.setVisible(true);

	}

	public void actionPerformed(ActionEvent actionEvent) {
		
		String yesCommand = new String("Yes");
		String noCommand = new String("No");
		String maybeCommand = new String("Maybe");		
		String backCommand = new String("Back");
		String openFilesCommand = new String("Load Files");
		String showLogCommand = new String("Show log");
		String quitCommand = new String("Quit");

		String command = actionEvent.getActionCommand();
		if (command.equals(noCommand)) {
			noButton.requestFocusInWindow();
			recordVerdict(-1);
			
		} else if (command.equals(maybeCommand)) {
			maybeButton.requestFocusInWindow();
			recordVerdict(0);
			
		} else if (command.equals(yesCommand)) {
			yesButton.requestFocusInWindow();
			recordVerdict(1);
			
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
			
			dcd.pack();
			dcd.setVisible(true);

			// check all the dialog fields are valid paths to data
			String bam = dcd.getBam().trim();
			ArrayList<String> messages = new ArrayList<String>();
			
			if (bam.equals("")) {
				messages.add("ERROR: BAM(s) file field empty");
			}
			
			String var = dcd.getVar();
			if (var.equals("")) {
				messages.add("ERROR: Variant file field empty");
			}
			
			String ref = dcd.getRef();
			if (ref.equals("")) {
				messages.add("ERROR: Reference file field empty");
			}
			
			if (messages.size() > 0) {
				// we have found an error so print it and return
				printMessage(messages);
				return;
			} else {				
				// if we fail to open any bams then go no further
				if (!openBams(bam)) {
					return;
				}
				// if we fail to load any variants then go no further
				if (!openVariants(var)) {
					return;
				}
				// if we fail to load a reference then go no further
				if (!openReference(ref)) {
					return;
				}
			}
			
			// if opening all the required files is successful then load the first variant
			currentVariantIndex = 0;
			displayVariantIndex = currentVariantIndex;
			activeScorePanel(true);
			drawReads();
				
		} else if (command.equals(quitCommand)) {
			System.exit(0);
		}
	}

	private boolean openReference(String ref) {
		try {
			File refFile = new File(ref);
			referenceFile = new IndexedFastaSequenceFile(refFile);
			refName = refFile.getName();
			printMessage("Loaded reference file: '"+ ref +"'");
		} catch (FileNotFoundException e) {
			JOptionPane.showMessageDialog(null, "Unable to open the selected reference file:\n" + e.toString(), "Error", JOptionPane.ERROR_MESSAGE);
			return false;
		}		
		return true;
	}

	private boolean openVariants(String var) {		
		VariantListFile variantListFile = null;
		try {
			variantListFile = new VariantListFile(var);
		} catch (IOException e) {
			JOptionPane.showMessageDialog(null, "Unable to open the selected variants file:\n" + e.toString(), "Error", JOptionPane.ERROR_MESSAGE);
			return false;
		}
		variantList = variantListFile.getVariants();
		printMessage(variantList.size() + " variants loaded from: '"+ var +"'");
		output = checkOverwriteFile(new File(var + ".scores"));
		
		if (output == null) {			
			JOptionPane.showMessageDialog(null, "Unable to open the selected output file:\n", "Error", JOptionPane.ERROR_MESSAGE);
			return false;
		}
		
		printMessage("Writing scores to: '"+ output.getAbsolutePath() +"'");
		return true;
	}

	private boolean openBams(String bam) {
		bams = new ArrayList<SAMFileReader>();
		bamNames = new ArrayList<String>();
		// if the bam path ends .bam just open the single file			
		if (bam.endsWith(".bam")) {
			File bamFile = new File(bam);			
			if (bamFile.exists()) {
				if (checkBamIndex(bam)) {
					bams.add(new SAMFileReader(bamFile));
					bamNames.add(bamFile.getName());
					printMessage("Loaded BAM: '"+ bam +"'");
				}
			} else {
				printMessage("Can't load BAM: '"+ bam +"'");
				return false;
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
						
						// check if the line contains 1 or 2 columns
						currentline = currentline.trim();
						String[] values = currentline.split("\\t+");
						String bamName = null;
						String bamDisplayName = null;
						
						if (values.length == 1) {
							bamName = values[0];
							bamDisplayName = values[0];
						} else if (values.length > 1) {
							// not expecting any more than 2 columns so ignore addition columns
							bamName = values[0];
							bamDisplayName = values[1];
						}
						
						String b = bamDir + File.separator + bamName;
						File bamFile = new File(b);
						if (bamFile.exists() && checkBamIndex(b)) {
							bams.add(new SAMFileReader(bamFile));
							bamNames.add(bamDisplayName);
							printMessage("Loaded BAM: '"+ b +"'");
						}
					}
				}
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		return true;
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
			
			JOptionPane.showMessageDialog(null, msg, "Error", JOptionPane.ERROR_MESSAGE);
			
			return Boolean.FALSE;
		} else {
			return Boolean.TRUE;	
		}
		
	}

	private void drawReads() {
		
		trackContainer.removeAll();
		nameContainer.removeAll();
		
		if (bams.size() > 0) {
			Variant v = variantList.get(displayVariantIndex);
			
			this.setTitle(v.getName());
			
			trackContainer.removeAll();
			nameContainer.removeAll();
			
			// calculate the height of each panel		
			int refPanelHeight = 50;
			int border = 40;
			int height = trackContainer.getHeight() - (controlsPanel.getHeight() + mb.getHeight()) - refPanelHeight - border;
			int trackHeight = height / bams.size();
			int width = trackContainer.getWidth();
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
			// add the variant base to the bases either side
			int refEnd = v.getPos() + basesOneSide + 1;

			ReferencePanel rp = new ReferencePanel(pixPerBase, v, referenceFile, refStart, refEnd);
			rp.setPreferredSize(new Dimension(width, refPanelHeight));
			trackContainer.add(rp);
						
			// add bam names
			for (String bamName : bamNames) {				
				
				JPanel nameHolder = new JPanel();
				nameHolder.add(new JLabel(bamName));
				nameHolder.setPreferredSize(new Dimension(100, trackHeight));
				nameHolder.setBorder(BorderFactory.createLineBorder(Color.black));
				nameContainer.add(nameHolder);
			}
			
			JPanel refNameHolder = new JPanel();
			refNameHolder.add(new JLabel(refName));
			refNameHolder.setPreferredSize(new Dimension(100, refPanelHeight));						
			nameContainer.add(refNameHolder);
			
			this.setVisible(true);
			this.repaint();
		}
		
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
			controlsPanel.setVisible(true);
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
			int n = JOptionPane.showConfirmDialog(	this.getContentPane(), "The file "	+ file.getName() + " already exists\n would you like to overwrite this file?", "Overwrite file?", JOptionPane.YES_NO_OPTION, JOptionPane.QUESTION_MESSAGE);
			// n 0 = yes 1 = no
			if (n == 1) {	
				JFileChooser jfc = new JFileChooser("user.dir");
				if (jfc.showSaveDialog(this) == JFileChooser.APPROVE_OPTION) {
					file = new File(jfc.getSelectedFile().getAbsolutePath());
				} else {
					// we have no output file so don't load the data
					file = null;
				}
			}
		}
		return file;
	}

	private void recordVerdict(int v) {
		if (isBackSNP()) {
			listScores.put(variantList.get(displayVariantIndex), v);
			try {
				printScores();
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			setBackSNP(false);
			if (displayVariantIndex == 0) {
				backButton.setEnabled(true);
			}
			displayVariantIndex = currentVariantIndex;
			drawReads();
		} else {
			listScores.put(variantList.get(currentVariantIndex), v);
			try {
				printScores();
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
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

	private void printMessage(String message) {
		
		SimpleDateFormat sdf = new SimpleDateFormat("h:mm a");
		String formattedDate = sdf.format(new Date());
		
		messageText.append(formattedDate + ": " + message + "\n");		
	}
	
	private void printMessage(ArrayList<String> messages) {
		
		SimpleDateFormat sdf = new SimpleDateFormat("h:mm a");
		String formattedDate = sdf.format(new Date());
		
		for (String message : messages) {
			messageText.append(formattedDate + ": " + message + "\n");
		}
		controlsPanel.repaint();
		controlsPanel.validate();
	}

	public void refreshPlot() {
		if (displayVariantIndex != null) {
			drawReads();
		}
	}

	@Override
	public void componentHidden(ComponentEvent e) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void componentMoved(ComponentEvent e) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void componentResized(ComponentEvent e) {
		if (displayVariantIndex != null) {
			drawReads();
		}		
	}

	@Override
	public void componentShown(ComponentEvent e) {
		// TODO Auto-generated method stub
		
	}

}
