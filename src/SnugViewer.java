import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JCheckBoxMenuItem;
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
import javax.swing.JSplitPane;
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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;

import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.ListIterator;
import java.util.ArrayList;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import com.jcraft.jsch.*;

public class SnugViewer extends JFrame implements ActionListener,
		ComponentListener {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1122484803356474269L;

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
	private JMenu viewMenu;
	private JCheckBoxMenuItem collapseViewCB;
	private DataConnectionDialog dcd;
	private JMenuBar mb;
	private JPanel trackContainer;
	private JPanel nameContainer;
	private String refName;
	private JPanel contentPanel;

	//
	// Global arrays
	//
	// for remote file access
	private Session remoteSession;
	// bamNames - array containing the samfilereader objects for each bam
	private ArrayList<SAMFileReader> bams;
	// bamNames - array containing the bam file name with full path removed uses
	// to display the bam name in the main display
	private ArrayList<String> bamNames;
	// remoteBams - array containing just the bam paths
	private ArrayList<String> remoteBams;

	//
	// Global variables
	//
	private boolean backSNP;
	private IndexedFastaSequenceFile referenceFile;
	// set the space given to each base as it is drawn
	private int pixPerBase = 15;
	// set the default view to be expanded
	private boolean collapseView = false;
	private int totalBams = 0;

	public static void main(String[] args) {
		new SnugViewer();
	}

	SnugViewer() {
		super("Snug...");

		dcd = new DataConnectionDialog(this);
		mb = new JMenuBar();

		int menumask = Toolkit.getDefaultToolkit().getMenuShortcutKeyMask();

		fileMenu = new JMenu("File");
		JMenuItem loadFiles = new JMenuItem("Load Files");
		loadFiles.setAccelerator(KeyStroke
				.getKeyStroke(KeyEvent.VK_B, menumask));
		loadFiles.addActionListener(this);
		fileMenu.add(loadFiles);

		mb.add(fileMenu);

		viewMenu = new JMenu("View");
		collapseViewCB = new JCheckBoxMenuItem("Collapse view");
		collapseViewCB.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_C,
				menumask));
		collapseViewCB.addActionListener(this);
		collapseViewCB.setEnabled(false);
		viewMenu.add(collapseViewCB);

		mb.add(viewMenu);

		if (!(System.getProperty("os.name").toLowerCase().contains("mac"))) {
			JMenuItem quitItem = new JMenuItem("Quit");
			quitItem.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_Q,
					menumask));
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
		scorePanel.registerKeyboardAction(this, "Yes",
				KeyStroke.getKeyStroke(KeyEvent.VK_Y, 0),
				JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT);
		scorePanel.add(noButton);
		scorePanel.registerKeyboardAction(this, "No",
				KeyStroke.getKeyStroke(KeyEvent.VK_N, 0),
				JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT);
		scorePanel.add(maybeButton);
		scorePanel.registerKeyboardAction(this, "Maybe",
				KeyStroke.getKeyStroke(KeyEvent.VK_M, 0),
				JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT);
		scorePanel.add(backButton);
		scorePanel.registerKeyboardAction(this, "Back",
				KeyStroke.getKeyStroke(KeyEvent.VK_B, 0),
				JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT);

		messageText = new JTextArea(4, 40);
		messageText.setEditable(false);
		messagePanel = new JScrollPane(messageText);

		controlsPanel = new JPanel();
		controlsPanel.add(scorePanel);
		controlsPanel.add(messagePanel);
		controlsPanel.setLayout(new BoxLayout(controlsPanel, BoxLayout.X_AXIS));
		controlsPanel.setMaximumSize(new Dimension(900, 200));
		controlsPanel.setVisible(false);
		controlsPanel.setAlignmentX(Component.LEFT_ALIGNMENT);

		trackContainer = new JPanel();
		trackContainer.setBackground(Color.white);
		trackContainer
				.setLayout(new BoxLayout(trackContainer, BoxLayout.Y_AXIS));
		trackContainer.setBorder(BorderFactory.createEmptyBorder(0, 0, 0, 10));

		nameContainer = new JPanel();

		nameContainer.setBackground(Color.white);
		nameContainer.setLayout(new BoxLayout(nameContainer, BoxLayout.Y_AXIS));
		nameContainer.setBorder(BorderFactory.createEmptyBorder(0, 10, 0, 0));

		JSplitPane mainPanel = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT,
				nameContainer, trackContainer);
		mainPanel.setOneTouchExpandable(true);
		mainPanel.setDividerLocation(100);
		mainPanel.setAlignmentX(Component.LEFT_ALIGNMENT);

		contentPanel = new JPanel();
		contentPanel.setPreferredSize(new Dimension(900, 600));
		contentPanel.setLayout(new BoxLayout(contentPanel, BoxLayout.Y_AXIS));
		contentPanel.add(mainPanel);
		contentPanel.add(controlsPanel);
		contentPanel.addComponentListener(this);

		this.setContentPane(contentPanel);
		this.pack();
		this.setVisible(true);

	}

	public void actionPerformed(ActionEvent actionEvent) {

		String yesCommand = "Yes";
		String noCommand = "No";
		String maybeCommand = "Maybe";
		String backCommand = "Back";
		String openFilesCommand = "Load Files";
		String quitCommand = "Quit";
		String collapseViewCommand = "Collapse view";

		String command = actionEvent.getActionCommand();
		if (command.equals(noCommand)) {
			noButton.requestFocusInWindow();
			try {
				recordVerdict(-1);
			} catch (JSchException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

		} else if (command.equals(maybeCommand)) {
			maybeButton.requestFocusInWindow();
			try {
				recordVerdict(0);
			} catch (JSchException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

		} else if (command.equals(yesCommand)) {
			yesButton.requestFocusInWindow();
			try {
				recordVerdict(1);
			} catch (JSchException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (IOException e) {
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
			try {
				drawReads();
			} catch (JSchException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		} else if (command.equals(collapseViewCommand)) {

			if (collapseViewCB.isSelected()) {
				collapseView = true;
			} else {
				collapseView = false;
			}
			if (displayVariantIndex != null) {
				try {
					drawReads();
				} catch (JSchException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		} else if (command.equals(openFilesCommand)) {
			dcd.pack();
			dcd.setVisible(true);

			if (dcd.isFilled()) {
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
					// TODO need to create an alert not print to the logging
					// panel
					// because the logging panel will not be visable
					// use logging panel for only success messages
					printMessage(messages);
					return;
				} else {

					remoteSession = null;
					// if any of the files have been selected to come from a
					// remote source then create a jsch session
					// otherwise the remoteSession will remain null
					if (dcd.remoteBam() | dcd.remoteRef() | dcd.remoteVar()) {
						try {
							remoteSession = createRemoteSession();
						} catch (JSchException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();

							// JOptionPane.showMessageDialog(
							// null,
							// "Unable to open the selected variants file:\n"
							// + e.toString(), "Error",
							// JOptionPane.ERROR_MESSAGE);

						}
					}

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

				// if opening all the required files is successful then load the
				// first variant
				currentVariantIndex = 0;
				displayVariantIndex = currentVariantIndex;
				activeScorePanel(true);
				try {
					drawReads();
				} catch (JSchException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}

		} else if (command.equals(quitCommand)) {
			System.exit(0);
		}
	}

	private Session createRemoteSession() throws JSchException {
		String host = dcd.getHost().trim();
		String user = dcd.getUsername().trim();
		String pass = dcd.getPassword();
		int port = dcd.getPort();

		java.util.Properties config = new java.util.Properties();
		config.put("StrictHostKeyChecking", "no");

		JSch jsch = new JSch();
		Session session = jsch.getSession(user, host, port);
		session.setConfig(config);
		session.setPassword(pass);
		session.connect();

		return session;
	}

	private boolean openReference(String ref) {
		try {
			File refFile = new File(ref);
			referenceFile = new IndexedFastaSequenceFile(refFile);
			refName = refFile.getName();
			printMessage("Loaded reference file: '" + ref + "'");
		} catch (FileNotFoundException e) {
			JOptionPane.showMessageDialog(
					null,
					"Unable to open the selected reference file:\n"
							+ e.toString(), "Error", JOptionPane.ERROR_MESSAGE);
			return false;
		}
		return true;
	}

	private boolean openVariants(String var) {
		VariantListFile variantListFile = null;
		try {
			variantListFile = new VariantListFile(var);
		} catch (IOException e) {
			JOptionPane.showMessageDialog(
					null,
					"Unable to open the selected variants file:\n"
							+ e.toString(), "Error", JOptionPane.ERROR_MESSAGE);
			return false;
		} catch (NumberFormatException e) {
			JOptionPane.showMessageDialog(null,
					"Position value in variants file is not a number", "Error",
					JOptionPane.ERROR_MESSAGE);
			return false;
		} catch (UnrecognisedVariantFileFormat e) {
			JOptionPane.showMessageDialog(null,
					"Unrecognised variant file format: " + e.getMessage(),
					"Error", JOptionPane.ERROR_MESSAGE);
			return false;
		}

		variantList = variantListFile.getVariants();
		printMessage(variantList.size() + " variants loaded from: '" + var
				+ "'");
		output = checkOverwriteFile(new File(var + ".scores"));

		if (output == null) {
			JOptionPane.showMessageDialog(null,
					"Unable to open the selected output file:\n", "Error",
					JOptionPane.ERROR_MESSAGE);
			return false;
		}

		printMessage("Writing scores to: '" + output.getAbsolutePath() + "'");
		return true;
	}

	private boolean openBams(String bam) {
		bams = new ArrayList<SAMFileReader>();
		remoteBams = new ArrayList<String>();
		bamNames = new ArrayList<String>();
		// if the bam path ends .bam just open the single file
		if (bam.endsWith(".bam")) {
			// if the bam file(s) are local the create an array on samfilereader
			// objects
			if (remoteSession == null) {
				File bamFile = new File(bam);
				if (bamFile.exists()) {
					if (checkBamIndex(bam)) {
						bams.add(new SAMFileReader(bamFile));
						bamNames.add(bamFile.getName());
						printMessage("Loaded BAM: '" + bam + "'");
					} else {
						return false;
					}
				} else {
					printMessage("Can't load BAM: '" + bam + "'");
					return false;
				}

				// if the bam file(s) are remote then just record the path to
				// each
			} else {
				// check for an index files
				// if an index can't be found then run samtools index
				remoteBams.add(bam);
				bamNames.add(bam);
				printMessage("Loaded BAM: '" + bam + "'");
			}
			totalBams = 1;

			// if the bam path ends .bams the file is assumed to contain a list
			// of bam files
		} else if (bam.endsWith(".bams")) {

			File bamFiles = new File(bam);
			String absolutePath = bamFiles.getAbsolutePath();
			String bamDir = absolutePath.substring(0, absolutePath.lastIndexOf(File.separator));

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
					// check that the line is not blank
					if (currentline == null || currentline.trim().equals("")) {
						continue;
					} else {

						// check if the line contains 1 or 2 columns
						// the second column can contain a short name for the bam
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

						if (remoteSession == null) {
							String b = bamDir + File.separator + bamName;
							File bamFile = new File(b);
							if (bamFile.exists() && checkBamIndex(b)) {
								bams.add(new SAMFileReader(bamFile));
								bamNames.add(bamDisplayName);
								printMessage("Loaded BAM: '" + b + "'");
							}
						} else {
							// assume to start with that the list of bams is local and
							// the paths in the file point to a remote location

							// therefore we dont care about the directory the file is in on
							// the local machine
							remoteBams.add(bamName);
							bamNames.add(bamDisplayName);
							printMessage("Loaded BAM: '" + bamName + "'");

						}						
					}
				}
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

			if (remoteSession == null) {
				totalBams = bams.size();
			} else {
				totalBams = remoteBams.size();
			}

		} else {
			JOptionPane.showMessageDialog(null,
					"Bam file should end .bam or .bams:\n", "Error",
					JOptionPane.ERROR_MESSAGE);
			return false;
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

			JOptionPane.showMessageDialog(null, msg, "Error",
					JOptionPane.ERROR_MESSAGE);

			return Boolean.FALSE;
		} else {
			return Boolean.TRUE;
		}

	}

	private void drawReads() throws JSchException, IOException {

		trackContainer.removeAll();
		nameContainer.removeAll();

		if (totalBams > 0) {
			Variant v = variantList.get(displayVariantIndex);

			this.setTitle(v.getName());

			trackContainer.removeAll();
			nameContainer.removeAll();

			// calculate the height of each panel
			int refPanelHeight = 50;
			int border = 40;
			int height = trackContainer.getHeight()
					- (controlsPanel.getHeight() + mb.getHeight())
					- refPanelHeight - border;
			int trackHeight = height / totalBams;
			int width = trackContainer.getWidth();
			int visableBases = width / pixPerBase;

			// if the track height is too low then add another column of
			// alignments

			// dont draw alignments if the height and width gets too low - warn
			// the user of too many bams?

			// minimum height and width values change depending if you are using
			// collapsed or expanded view

			// assume the middle base is the variant - how many bases are either
			// side of it
			int basesOneSide = (visableBases - 1) / 2;

			// distance in pix to the variant site
			int offset = basesOneSide * pixPerBase;

			// get the reads for the region
			ArrayList<SAMRecord> reads = new ArrayList<SAMRecord>();

			if (remoteSession == null) {
				for (SAMFileReader sfr : bams) {

					SAMRecordIterator readsItr = sfr.queryOverlapping(
							v.getChr(), v.getPos(), v.getPos());

					while (readsItr.hasNext()) {
						reads.add(readsItr.next());
					}
					readsItr.close();

					AlignmentPanel ap = new AlignmentPanel(pixPerBase, offset,
							reads, v, height, referenceFile, collapseView);

					ap.setPreferredSize(new Dimension(width, trackHeight));
					setBackground(Color.white);
					ap.setBorder(BorderFactory.createLineBorder(Color.black));

					trackContainer.add(ap);
				}
			} else {
				for (String bam : remoteBams) {
					String cmd = "samtools view -h " + bam + " " + v.getChr()
							+ ":" + v.getPos() + "-" + v.getPos();
					// create a channel connected to a remotely executing
					// program
					Channel channel = remoteSession.openChannel("exec");

					((ChannelExec) channel).setCommand(cmd);

					// capture the error stream and parse it to deliver sensible
					// message to the user
					// The standard error output of the remote process will be
					// sent to this stream
					((ChannelExec) channel).setErrStream(System.err);

					InputStream inputstream = channel.getInputStream();

					channel.connect();

					if (channel.getExitStatus() != 0) {
						SAMRecordIterator readsItr = new SAMFileReader(
								inputstream, true).iterator();

						while (readsItr.hasNext()) {
							reads.add(readsItr.next());
						}
					} else {
						// tell the user that the remote command failed
					}
					
					AlignmentPanel ap = new AlignmentPanel(pixPerBase, offset,
							reads, v, height, referenceFile, collapseView);

					ap.setPreferredSize(new Dimension(width, trackHeight));
					setBackground(Color.white);
					ap.setBorder(BorderFactory.createLineBorder(Color.black));

					trackContainer.add(ap);
					
				}
				
				
				
				
			}

			AlignmentPanel ap = new AlignmentPanel(pixPerBase, offset, reads,
					v, height, referenceFile, collapseView);

			ap.setPreferredSize(new Dimension(width, trackHeight));
			setBackground(Color.white);
			ap.setBorder(BorderFactory.createLineBorder(Color.black));

			trackContainer.add(ap);

			// draw the reference panel
			int refStart = v.getPos() - basesOneSide;
			// add the variant base to the bases either side
			int refEnd = v.getPos() + basesOneSide + 1;

			ReferencePanel rp = new ReferencePanel(pixPerBase, v,
					referenceFile, refStart, refEnd);
			rp.setPreferredSize(new Dimension(width, refPanelHeight));
			trackContainer.add(rp);

			// add bam names
			for (String bamName : bamNames) {
				JPanel nameHolder = new JPanel();
				nameHolder.add(new JLabel(bamName));
				nameHolder.setPreferredSize(new Dimension(100, trackHeight));
				nameHolder.setBorder(BorderFactory
						.createLineBorder(Color.black));
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
			collapseViewCB.setEnabled(true);
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
			collapseViewCB.setEnabled(false);
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

	private void recordVerdict(int v) throws JSchException, IOException {
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

	public void refreshPlot() throws JSchException, IOException {
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
			try {
				drawReads();
			} catch (JSchException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			} catch (IOException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
		}
	}

	@Override
	public void componentShown(ComponentEvent e) {
		// TODO Auto-generated method stub

	}

}
