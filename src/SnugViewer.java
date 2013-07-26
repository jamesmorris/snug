import javax.imageio.ImageIO;
import javax.swing.*;
import javax.swing.border.Border;
import javax.swing.border.EmptyBorder;
import javax.swing.border.LineBorder;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.picard.reference.IndexedFastaSequenceFile;

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

public class SnugViewer extends JFrame implements ActionListener {

	String plottedSNP = null;
	private File output;
	private ArrayList<Variant> variantList;
	private Variant currentVariantinList;
	private Variant currentVariant = null;
	private int currentVariantIndex;
	private int displayVariantIndex;
	private boolean backVariant = false;
	private boolean bamFile = false;
	private boolean variantFile = false;
	private boolean endOfList;
	private HashMap<Variant, Integer> listScores = new HashMap<Variant, Integer>();
	JFileChooser jfc;
	private JPanel scorePanel;
	private JPanel messagePanel;
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
	
	JPanel contentPanel;

	SAMFileReader bamReader;

	public enum MouseMode {
		LASSO, ZOOM
	};

	private MouseMode mouseMode;
	private boolean backSNP;
	private AlignmentPanel ap;
	private IndexedFastaSequenceFile referenceFile;

	public static void main(String[] args) {

		new SnugViewer();

	}

	SnugViewer() {
		super("Snug...");

		dcd = new DataConnectionDialog(this);
		JMenuBar mb = new JMenuBar();

		int menumask = Toolkit.getDefaultToolkit().getMenuShortcutKeyMask();

		fileMenu = new JMenu("File");
		loadFiles = new JMenuItem("Load Files");
		loadFiles.setAccelerator(KeyStroke
				.getKeyStroke(KeyEvent.VK_B, menumask));
		loadFiles.addActionListener(this);
		fileMenu.add(loadFiles);

		mb.add(fileMenu);

		if (!(System.getProperty("os.name").toLowerCase().contains("mac"))) {
			JMenuItem quitItem = new JMenuItem("Quit");
			quitItem.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_Q,
					menumask));
			quitItem.addActionListener(this);
			fileMenu.add(quitItem);
		}

		JMenu logMenu = new JMenu("Log");
		showLogItem = new JMenuItem("Show log");
		showLogItem.addActionListener(this);
		logMenu.add(showLogItem);
		mb.add(logMenu);

		setJMenuBar(mb);

		JPanel controlsPanel = new JPanel();

		scorePanel = new JPanel();
		scorePanel.add(new JLabel("Approve?"));

		yesButton = new JButton("Yes");
		scorePanel.registerKeyboardAction(this, "Yes",
				KeyStroke.getKeyStroke(KeyEvent.VK_Y, 0),
				JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT);
		yesButton.addActionListener(this);
		yesButton.setEnabled(false);
		scorePanel.add(yesButton);

		maybeButton = new JButton("Maybe");
		scorePanel.registerKeyboardAction(this, "Maybe",
				KeyStroke.getKeyStroke(KeyEvent.VK_M, 0),
				JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT);
		maybeButton.addActionListener(this);
		maybeButton.setEnabled(false);
		scorePanel.add(maybeButton);

		noButton = new JButton("No");
		scorePanel.registerKeyboardAction(this, "No",
				KeyStroke.getKeyStroke(KeyEvent.VK_N, 0),
				JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT);
		noButton.addActionListener(this);
		noButton.setEnabled(false);
		scorePanel.add(noButton);

		backButton = new JButton("Back");
		scorePanel.registerKeyboardAction(this, "Back",
				KeyStroke.getKeyStroke(KeyEvent.VK_B, 0),
				JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT);
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

		controlsPanel.setMaximumSize(new Dimension(2000, (int) controlsPanel
				.getPreferredSize().getHeight()));
		controlsPanel.setMinimumSize(new Dimension(10, (int) controlsPanel
				.getPreferredSize().getHeight()));

		mouseMode = MouseMode.ZOOM;

		ap = new AlignmentPanel();
		
		contentPanel = new JPanel();
		contentPanel.setLayout(new BoxLayout(contentPanel, BoxLayout.Y_AXIS));
		contentPanel.add(ap);
		contentPanel.add(controlsPanel);

		addWindowListener(new WindowAdapter() {

			public void windowClosing(WindowEvent e) {
				System.exit(0);
			}
		});

		this.setContentPane(contentPanel);
		this.pack();
		this.setVisible(true);

		// add a panel and to it add a text string for test
		ap.repaint();

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

			dcd.pack();
			dcd.setVisible(true);

			// check all the dialog fields are valid paths to data
			
			String bam = dcd.getBam();
			File bamFile = new File(bam);

			if (bamFile.exists()) {
				bamReader = new SAMFileReader(bamFile);
			}

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
			} else {
				ld.log("BAM loaded");
			}

			currentVariantIndex = 0;
			displayVariantIndex = currentVariantIndex;
			String var = dcd.getVar();
			VariantListFile variantListFile = new VariantListFile(var);
			variantList = variantListFile.getVariants();
			ld.log(variantList.size() + " variants loaded");
			output = checkOverwriteFile(new File(var + ".scores"));
			
			try {
				ap.setReference(new IndexedFastaSequenceFile(new File(dcd.getRef())));
				ld.log("Reference file loaded");
			} catch (FileNotFoundException e) {
				ld.log("Reference file not found");
			}
			
			yesButton.setEnabled(true);
			noButton.setEnabled(true);
			maybeButton.setEnabled(true);
			yesButton.requestFocusInWindow();
			drawReads();			
				
		} else if (command.equals(quitCommand)) {
			System.exit(0);
		} else if (command.equals(showLogCommand)) {
			ld.setVisible(true);
		}
	}

	private void drawReads() {
		Variant v = variantList.get(displayVariantIndex);
		ArrayList<SAMRecord> reads = new ArrayList<SAMRecord>();		
		SAMRecordIterator readsItr = bamReader.queryOverlapping(v.getChr(), v.getPos(), v.getPos());
		
        while (readsItr.hasNext()) {	
        	reads.add(readsItr.next());  	
        }
        readsItr.close();
        
        Integer region_start = null;
		Integer region_end = 0;
		
		for (SAMRecord read : reads) {
			if (region_start == null) {
				region_start = read.getAlignmentStart();
			}
			if (read.getAlignmentEnd() > region_end) {
				region_end = read.getAlignmentEnd();
			}
		}
		
		ap.setReads(reads);
		ap.setVariant(v);
		ap.repaint();

		this.pack();
		this.setVisible(true);
		this.setTitle(v.getName());

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

	protected MouseMode getMouseMode() {
		return mouseMode;
	}

	protected void setMouseMode(MouseMode mm) {
		this.mouseMode = mm;
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
		if (currentVariant != null) {
			drawReads();
		}
	}
}
