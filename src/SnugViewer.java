import javax.swing.BorderFactory;
import javax.swing.BoundedRangeModel;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
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
import javax.swing.JRadioButtonMenuItem;
import javax.swing.JScrollBar;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTextArea;
import javax.swing.KeyStroke;
import javax.swing.ScrollPaneConstants;
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
import java.io.OutputStream;
import java.io.PrintStream;

import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.HashMap;
import java.util.ListIterator;
import java.util.ArrayList;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFormatException;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import com.jcraft.jsch.*;

public class SnugViewer extends JFrame implements ActionListener, ComponentListener {

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
	// set the height of the reference sequence panel
	int refPanelHeight = 60;
	// set the default view settings
	private boolean collapseView = false;
	private boolean refernceOnTop = false;
	private boolean hideRefBases = false;
	private boolean mappingQuality = false;
	private int totalBams = 0;
	
	//
	// UI commands
	//
	private String yesCommand = "Yes";
	private String noCommand = "No";
	private String maybeCommand = "Maybe";
	private String backCommand = "Back";
	private String openFilesCommand = "Load Files";
	private String quitCommand = "Quit";
	private String mappingQualityCommand = "Colour by mapping quality";
	private String collapseAllCommand = "All";
	private String collapseRefCommand = "Reference";
	private String collapseNoneCommand = "None";
	private String refTopCommand = "Top";
	private String refBotCommand = "Bottom";
	private String colourBaseCommand = "By Base Qual";
	private String colourReadCommand = "By Read Qual";

	private JMenu collapseSubmenu;
	private JMenu referenceSubMenu;
	private JMenu colourSubMenu;
	private JRadioButtonMenuItem collapseAllButton;
	private JRadioButtonMenuItem collapseRefButton;
	private JRadioButtonMenuItem collapseNoneButton;
	private JRadioButtonMenuItem refTopButton;
	private JRadioButtonMenuItem refBotButton;
	private JRadioButtonMenuItem colourBaseButton;
	private JRadioButtonMenuItem colourReadButton;
	
	
	public static void main(String[] args) {
		new SnugViewer();
	}

	SnugViewer() {
		super("Snug...");

		dcd = new DataConnectionDialog(this);
		mb = new JMenuBar();

		int menumask = Toolkit.getDefaultToolkit().getMenuShortcutKeyMask();

		fileMenu = new JMenu("File");
		JMenuItem loadFiles = new JMenuItem(openFilesCommand);
		loadFiles.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_B, menumask));
		loadFiles.addActionListener(this);
		fileMenu.add(loadFiles);

		mb.add(fileMenu);

		viewMenu = new JMenu("View");

		collapseSubmenu = new JMenu("Collapse");
		
		ButtonGroup collapseGroup = new ButtonGroup();
		
		collapseAllButton = new JRadioButtonMenuItem(collapseAllCommand);
		collapseAllButton.addActionListener(this);
		collapseGroup.add(collapseAllButton);
		collapseSubmenu.add(collapseAllButton);
		
		collapseRefButton = new JRadioButtonMenuItem(collapseRefCommand);
		collapseRefButton.addActionListener(this);
		collapseGroup.add(collapseRefButton);
		collapseSubmenu.add(collapseRefButton);
		
		collapseNoneButton = new JRadioButtonMenuItem(collapseNoneCommand);
		collapseNoneButton.addActionListener(this);
		collapseGroup.add(collapseNoneButton);
		collapseSubmenu.add(collapseNoneButton);
		
		viewMenu.add(collapseSubmenu);		
		viewMenu.addSeparator();
		
		referenceSubMenu = new JMenu("Reference");
		
		ButtonGroup referenceGroup = new ButtonGroup();
		
		refTopButton = new JRadioButtonMenuItem(refTopCommand);
		refTopButton.addActionListener(this);
		referenceGroup.add(refTopButton);
		referenceSubMenu.add(refTopButton);
		
		refBotButton = new JRadioButtonMenuItem(refBotCommand);
		refBotButton.addActionListener(this);
		referenceGroup.add(refBotButton);
		referenceSubMenu.add(refBotButton);
		
		viewMenu.add(referenceSubMenu);
		viewMenu.addSeparator();

		colourSubMenu = new JMenu("Colour");
		
		ButtonGroup colourGroup = new ButtonGroup();
		
		colourBaseButton = new JRadioButtonMenuItem(colourBaseCommand);
		colourBaseButton.addActionListener(this);
		colourGroup.add(colourBaseButton);
		colourSubMenu.add(colourBaseButton);
		
		colourReadButton = new JRadioButtonMenuItem(colourReadCommand);
		colourReadButton.addActionListener(this);
		colourGroup.add(colourReadButton);
		colourSubMenu.add(colourReadButton);
		
		viewMenu.add(colourSubMenu);
		viewMenu.addSeparator();
			
		mb.add(viewMenu);

		if (!(System.getProperty("os.name").toLowerCase().contains("mac"))) {
			JMenuItem quitItem = new JMenuItem("Quit");
			quitItem.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_Q, menumask));
			quitItem.addActionListener(this);
			fileMenu.add(quitItem);
		}

		setJMenuBar(mb);

		yesButton = new JButton(yesCommand);
		yesButton.addActionListener(this);
		yesButton.setEnabled(false);

		maybeButton = new JButton(maybeCommand);
		maybeButton.addActionListener(this);
		maybeButton.setEnabled(false);

		noButton = new JButton(noCommand);
		noButton.addActionListener(this);
		noButton.setEnabled(false);

		backButton = new JButton(backCommand);
		backButton.addActionListener(this);
		backButton.setEnabled(false);

		scorePanel = new JPanel();
		scorePanel.add(new JLabel("Approve?"));
		scorePanel.add(yesButton);
		scorePanel.registerKeyboardAction(this, yesCommand, KeyStroke.getKeyStroke(KeyEvent.VK_Y, 0), JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT);
		scorePanel.add(noButton);
		scorePanel.registerKeyboardAction(this, noCommand, KeyStroke.getKeyStroke(KeyEvent.VK_N, 0), JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT);
		scorePanel.add(maybeButton);
		scorePanel.registerKeyboardAction(this, maybeCommand, KeyStroke.getKeyStroke(KeyEvent.VK_M, 0), JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT);
		scorePanel.add(backButton);
		scorePanel.registerKeyboardAction(this, backCommand, KeyStroke.getKeyStroke(KeyEvent.VK_B, 0), JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT);

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
		trackContainer.setLayout(new BoxLayout(trackContainer, BoxLayout.Y_AXIS));
		trackContainer.setBorder(BorderFactory.createEmptyBorder(0, 0, 10, 10));

		nameContainer = new JPanel();
		nameContainer.setBackground(Color.white);
		nameContainer.setLayout(new BoxLayout(nameContainer, BoxLayout.Y_AXIS));
		nameContainer.setBorder(BorderFactory.createEmptyBorder(0, 10, 10, 0));

		JSplitPane mainPanel = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, nameContainer, trackContainer);
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
		} else if (command.equals(collapseAllCommand)) {
			collapseView = true;
			hideRefBases = false;
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
		} else if (command.equals(collapseRefCommand)) {
			hideRefBases = true;
			collapseView = false;
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
			
		} else if (command.equals(collapseNoneCommand)) {
			hideRefBases = false;
			collapseView = false;
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
			
		} else if (command.equals(refTopCommand)) {
			refernceOnTop = true;
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
			
		} else if (command.equals(refBotCommand)) {
			refernceOnTop = false;
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
			
		} else if (command.equals(colourBaseCommand)) {
			mappingQuality = false;
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
			
		} else if (command.equals(colourReadCommand)) {
			mappingQuality = true;
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
			
		}
		
		else if (command.equals(openFilesCommand)) {
			dcd.pack();
			dcd.setVisible(true);

			if (dcd.isFilled()) {
				
				boolean bamSuccess;
				boolean varSuccess;
				boolean refSuccess;
				// Initialise the arrays that will hold the bam file paths
				bams = new ArrayList<SAMFileReader>();
				remoteBams = new ArrayList<String>();
				bamNames = new ArrayList<String>();

				remoteSession = null;
				// if any of the files have been selected to come from a remote
				// source then create a jsch session
				if (dcd.remoteBam() | dcd.remoteRef() | dcd.remoteVar()) {
					if (dcd.getHost().equals("")) {
						JOptionPane.showMessageDialog(null, "Loading files failed: Host name missing, please try again", "Error", JOptionPane.ERROR_MESSAGE);
						return;
					} else if (dcd.getUsername().equals("")) {
						JOptionPane.showMessageDialog(null, "Loading files failed: User name missing, please try again", "Error", JOptionPane.ERROR_MESSAGE);
						return;
					} else if (dcd.getPassword().equals("")) {
						JOptionPane.showMessageDialog(null, "Loading files failed: Password missing, please try again", "Error", JOptionPane.ERROR_MESSAGE);
						return;
					} else {
						// if all the ssh details are filled out then try and
						// create a jsch session object
						try {
							remoteSession = createRemoteSession();
						} catch (JSchException e) {
							JOptionPane.showMessageDialog(null, "Connection to remote server failed:\n" + e.toString() + "\nPlease check your login details and try again", "Error", JOptionPane.ERROR_MESSAGE);
							return;
						}
					}
				}

				if (dcd.singleBam()) {
					if (dcd.getBam().equals("")) {
						JOptionPane.showMessageDialog(null, "Loading files failed: BAM file field empty, please try again", "Error", JOptionPane.ERROR_MESSAGE);
						return;
					} else {
						bamSuccess = openBam(dcd.getBam());
					}
				} else {
					if (dcd.getBamList().equals("")) {
						JOptionPane.showMessageDialog(null, "Loading files failed: BAM list field empty, please try again", "Error", JOptionPane.ERROR_MESSAGE);
						return;
					} else {
						bamSuccess = openBamList(dcd.getBamList());
					}
				}

				if (dcd.getVar().equals("")) {
					JOptionPane.showMessageDialog(null, "Loading files failed: Variant file field empty, please try again", "Error", JOptionPane.ERROR_MESSAGE);
					return;
				} else {
					varSuccess = openVariants(dcd.getVar());
				}

				if (dcd.getRef().equals("")) {
					JOptionPane.showMessageDialog(null, "Loading files failed: Reference file field empty, please try again", "Error", JOptionPane.ERROR_MESSAGE);
					return;
				} else {
					refSuccess = openReference(dcd.getRef());
				}

				if (bamSuccess && varSuccess && refSuccess) {
					// if opening all the required files is successful then load the first variant
					currentVariantIndex = 0;
					displayVariantIndex = currentVariantIndex;
					activeScorePanel(true);
					try {
						drawReads();
					} catch (JSchException e) {
						JOptionPane.showMessageDialog(null, e.toString(), "Error", JOptionPane.ERROR_MESSAGE);
					} catch (IOException e) {
						JOptionPane.showMessageDialog(null, e.toString(), "Error", JOptionPane.ERROR_MESSAGE);
					}
				} else {
					if (!bamSuccess) {
						JOptionPane.showMessageDialog(null, "Unable to open the bam file", "Error", JOptionPane.ERROR_MESSAGE);
					} else if (!varSuccess) {
						JOptionPane.showMessageDialog(null, "Unable to open the variants file", "Error", JOptionPane.ERROR_MESSAGE);
					} else {
						JOptionPane.showMessageDialog(null, "Unable to open the reference file", "Error", JOptionPane.ERROR_MESSAGE);
					}
				}
			}

		} else if (command.equals(quitCommand)) {
			System.exit(0);
		}
	}

	private boolean openBamList(String bamList) {

		File bamFiles = new File(bamList);
		String absolutePath = bamFiles.getAbsolutePath();
		String bamDir = absolutePath.substring(0, absolutePath.lastIndexOf(File.separator));

		// open the bams file and loop through all the bam files
		BufferedReader listReader = null;
		try {
			listReader = new BufferedReader(new FileReader(bamFiles));
		} catch (FileNotFoundException e) {
			JOptionPane.showMessageDialog(null, "Unable to open the bam list: " + e.toString(), "Error", JOptionPane.ERROR_MESSAGE);
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
					String[] values = currentline.trim().split("\\t+");
					String bamName = null;
					String bamDisplayName = null;

					if (values.length == 1) {
						bamName = values[0];
						bamDisplayName = values[0];
					} else if (values.length > 1) {
						// not expecting any more than 2 columns so ignore any
						// addition columns
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
						// assume to start with that the list of bams is local
						// and the paths in the file point to a remote location
						// therefore we dont care about the directory the file
						// is in on the local machine
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
		return true;
	}

	private boolean openBam(String bam) {

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
			// if the bam file is remote then just record the remote path to the
			// file
		} else {
			// TODO check for an index files if an index can't be found then run
			// samtools index
			remoteBams.add(bam);
			bamNames.add(bam);
			printMessage("Loaded BAM: '" + bam + "'");
		}
		totalBams = 1;

		return true;
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
		} catch (NumberFormatException e) {
			JOptionPane.showMessageDialog(null, "Position value in variants file is not a number", "Error", JOptionPane.ERROR_MESSAGE);
			return false;
		} catch (UnrecognisedVariantFileFormat e) {
			JOptionPane.showMessageDialog(null, "Unrecognised variant file format: " + e.getMessage(), "Error", JOptionPane.ERROR_MESSAGE);
			return false;
		}

		variantList = variantListFile.getVariants();
		printMessage(variantList.size() + " variants loaded from: '" + var + "'");
		output = checkOverwriteFile(new File(var + ".scores"));

		if (output == null) {
			JOptionPane.showMessageDialog(null, "Unable to open the selected output file:\n", "Error", JOptionPane.ERROR_MESSAGE);
			return false;
		}

		printMessage("Writing scores to: '" + output.getAbsolutePath() + "'");
		return true;
	}

	private Boolean checkBamIndex(String bam) {

		File bamIndexFile = new File(bam + ".bai");

		if (!bamIndexFile.exists()) {
			String ls = System.getProperty("line.separator");
			String msg = "BAM index file is missing. The BAM file needs to be sorted and indexed" + ls + "This can be done using samtools (http://samtools.sf.net/):" + ls + ls + "samtools sort <in.bam> <out.prefix>" + ls + "samtools index <sorted.bam>";

			JOptionPane.showMessageDialog(null, msg, "Error", JOptionPane.ERROR_MESSAGE);

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

			// calculate the height of each panel
			int trackHeight = (trackContainer.getHeight() - (controlsPanel.getHeight() + mb.getHeight()) - refPanelHeight) / totalBams;
			
			int alignmentStart = 0;
			int alignmentEnd = 0;
			boolean firstRead = true;
			
			// first loop through all the reads to get the size of the alignment
			ArrayList<ArrayList<SAMRecord>> readArray = new ArrayList<ArrayList<SAMRecord>>();
			
			if (remoteSession == null) {
				
				for (SAMFileReader sfr : bams) {
					ArrayList<SAMRecord> reads = new ArrayList<SAMRecord>();
					
					SAMRecordIterator readsItr = null;
					try {
						readsItr = sfr.queryOverlapping(v.getChr(), v.getPos(), v.getPos());
					} catch (SAMFormatException sfe) {
						JOptionPane.showMessageDialog(null, sfe.getMessage(), "Error", JOptionPane.ERROR_MESSAGE);
					}
										
					if (readsItr != null) {
						while (readsItr.hasNext()) {
							SAMRecord read = readsItr.next();
							if (firstRead) {
								alignmentStart = read.getAlignmentStart();
								alignmentEnd = read.getAlignmentEnd();
								firstRead = false;
							} else {
								if (read.getAlignmentStart() < alignmentStart) {
									alignmentStart = read.getAlignmentStart();
								}
								if (read.getAlignmentEnd() > alignmentEnd) {
									alignmentEnd = read.getAlignmentEnd();
								}
							}
							reads.add(read);
						}
						readsItr.close();
						readArray.add(reads);
					}					
				}
			} else {
				for (String bam : remoteBams) {
					String cmd = "samtools view -h " + bam + " " + v.getChr() + ":" + v.getPos() + "-" + v.getPos();
					// create a channel connected to a remotely executing program
					Channel channel = remoteSession.openChannel("exec");
					((ChannelExec) channel).setCommand(cmd);
					// The standard error output of the remote process will be sent to this stream
					OutputStream stderr = null;
					((ChannelExec) channel).setErrStream(stderr);
					
					InputStream inputstream = channel.getInputStream();
					InputStream errstream = ((ChannelExec) channel).getErrStream();
					channel.connect();
					
					ArrayList<SAMRecord> reads = new ArrayList<SAMRecord>();

					if (channel.getExitStatus() != 0) {
						SAMRecordIterator readsItr = new SAMFileReader(inputstream, true).iterator();
						while (readsItr.hasNext()) {
							SAMRecord read = readsItr.next();
							if (firstRead) {
								alignmentStart = read.getAlignmentStart();
								alignmentEnd = read.getAlignmentEnd();
								firstRead = false;
							} else {
								if (read.getAlignmentStart() < alignmentStart) {
									alignmentStart = read.getAlignmentStart();
								}
								if (read.getAlignmentEnd() > alignmentEnd) {
									alignmentEnd = read.getAlignmentEnd();
								}
							}
							reads.add(read);							
						}
					} else {
						JOptionPane.showMessageDialog(null, "Remote command failed: " + cmd, "Error", JOptionPane.ERROR_MESSAGE);
					}
					
					channel.disconnect();
					
					// convert the inputstream containing the standard error output into a string
					java.util.Scanner errString = new java.util.Scanner(errstream).useDelimiter("\\A");
				    
					if (errString.hasNext()) {
						JOptionPane.showMessageDialog(null, "Remote command failed: " + cmd + "\n" +errString.next(), "Error", JOptionPane.ERROR_MESSAGE);						
					}
					readArray.add(reads);					
				}
			}
		
			// calculate the size of the alignment
			int alignmentWidth = (alignmentEnd - alignmentStart) * pixPerBase;
			// distance in pix to the variant site from the start of the alignment
			int offset = (v.getPos() - alignmentStart) * pixPerBase;

			ReferencePanel rp = new ReferencePanel(pixPerBase, v, referenceFile, alignmentStart, alignmentEnd);
			rp.setPreferredSize(new Dimension(alignmentWidth, refPanelHeight));
			JScrollPane refScrollPane = new JScrollPane(rp);
			refScrollPane.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);
			JScrollBar refScrollBar = refScrollPane.getHorizontalScrollBar();
			BoundedRangeModel refScrollModel = refScrollBar.getModel();
			
			if (refernceOnTop) {
				trackContainer.add(refScrollPane);
			}
			
			// loop through each read array and generate the alignment panels
			for (ArrayList<SAMRecord> reads : readArray) {
				generateAlignmentPanel(pixPerBase, offset, reads, v, trackHeight, alignmentWidth, refScrollModel);
			}
			
			// if ref top is false then print the reference panel here
			if (!refernceOnTop) {
				trackContainer.add(refScrollPane);
			}
			
			// create a panel to hold the scrollbar
			JPanel sp = new JPanel();
			sp.setPreferredSize(new Dimension(alignmentWidth, 0));
			JScrollPane spp = new JScrollPane(sp);
			JScrollBar scrollBar = spp.getHorizontalScrollBar();
			scrollBar.setModel(refScrollModel);
			trackContainer.add(spp);
			
			// add track names
			JPanel refNameHolder = new JPanel();
			refNameHolder.add(new JLabel(refName));
			refNameHolder.setPreferredSize(new Dimension(100, refPanelHeight));
			
			if (refernceOnTop) {				
				nameContainer.add(refNameHolder);
			}
			
			for (String bamName : bamNames) {
				JPanel nameHolder = new JPanel();
				nameHolder.add(new JLabel(bamName));
				nameHolder.setPreferredSize(new Dimension(100, trackHeight));
				nameHolder.setBorder(BorderFactory.createLineBorder(Color.black));
				nameContainer.add(nameHolder);
			}
			
			if (!refernceOnTop) {
				nameContainer.add(refNameHolder);
			}
			
			JPanel nameScrollPanel = new JPanel();
			nameScrollPanel.setPreferredSize(new Dimension(100, 0));
			JScrollPane nameScrollPane = new JScrollPane(nameScrollPanel);
			nameScrollPane.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);
			nameScrollPane.setVerticalScrollBarPolicy(ScrollPaneConstants.VERTICAL_SCROLLBAR_NEVER);
			nameContainer.add(nameScrollPane);

			this.setVisible(true);
			
			// now everything is painted move the scrollbar to the middle of the window
			int mid = ((scrollBar.getMaximum() - scrollBar.getVisibleAmount()) - scrollBar.getMinimum()) / 2;
			scrollBar.setValue(mid);
			
			this.repaint();
		}

	}

	private void generateAlignmentPanel(int pixPerBase2, int offset, ArrayList<SAMRecord> reads, Variant v, int height, int width, BoundedRangeModel model) {
		if (reads.size() > 0) {
			AlignmentPanel ap = new AlignmentPanel(pixPerBase, offset, reads, v, referenceFile, collapseView, hideRefBases, mappingQuality);
			ap.setPreferredSize(new Dimension(width, height));
			setBackground(Color.white);
			ap.setBorder(BorderFactory.createLineBorder(Color.black));
			JScrollPane scrollPane = new JScrollPane(ap);
			scrollPane.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);
			scrollPane.setVerticalScrollBarPolicy(ScrollPaneConstants.VERTICAL_SCROLLBAR_AS_NEEDED);
			JScrollBar scrollBar = scrollPane.getHorizontalScrollBar();
			scrollBar.setModel(model);
			trackContainer.add(scrollPane);  			
		} else {
			JPanel jp = new JPanel();
			jp.setPreferredSize(new Dimension(width, height));
			jp.add(new JLabel("No reads found"));
			setBackground(Color.white);
			jp.setBorder(BorderFactory.createLineBorder(Color.black));
			trackContainer.add(jp);
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
			viewMenu.setEnabled(true);
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
			viewMenu.setEnabled(false);
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
			int n = JOptionPane.showConfirmDialog(this.getContentPane(), "The file " + file.getName() + " already exists\n would you like to overwrite this file?", "Overwrite file?", JOptionPane.YES_NO_OPTION, JOptionPane.QUESTION_MESSAGE);
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
				int n = JOptionPane.showConfirmDialog(this.getContentPane(), "Would you like to finish scoring the current list?", "Finish scoring list?", JOptionPane.YES_NO_OPTION, JOptionPane.QUESTION_MESSAGE);
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

	public void refreshPlot() throws JSchException, IOException {
		if (displayVariantIndex != null) {
			drawReads();
		}
	}

	@Override
	public void componentHidden(ComponentEvent e) {
		// required method override
	}

	@Override
	public void componentMoved(ComponentEvent e) {
		// required method override
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
		// required method override
	}

}
