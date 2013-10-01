import javax.swing.*;

import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;

public class DataConnectionDialog extends JDialog implements ActionListener {
    private String bamFile;
    private String bamList;
    private String varFile;
    private String refFile;
    private String host;
    private String port;
    private String username;
    private char[] password;
    private JTextField bamField;
    private JTextField bamListField;
    private JTextField varField;
    private JTextField refField;
    private JTextField hostField;
	private JTextField portField;
	private JTextField userField;
	private JPasswordField pf;
	JRadioButton bamLocalButton;
	JRadioButton bamRemoteButton;
	JRadioButton bamListLocalButton;
	JRadioButton bamListRemoteButton;
	JRadioButton varLocalButton;
	JRadioButton varRemoteButton;
	JRadioButton refLocalButton;
	JRadioButton refRemoteButton;
	JRadioButton singleButton;
    JRadioButton multipleButton;
	private boolean filled;
        
    JButton bamButton;
    JButton bamListButton;
    JButton varButton;
    JButton refButton;

    public DataConnectionDialog(JFrame parent){
        super(parent,"Data Connection",true);
        filled = false;
        JPanel contents = new JPanel();
        contents.setLayout(new BoxLayout(contents,BoxLayout.Y_AXIS));
        
        JPanel bamPanel = new JPanel();
        singleButton = new JRadioButton("");
        singleButton.setActionCommand("singleBam");
        singleButton.addActionListener(this);
        singleButton.setSelected(true);
        bamPanel.add(singleButton);
        bamPanel.add(new JLabel("BAM file: "));
        bamField = new JTextField(20);
        bamPanel.add(bamField);
        bamButton = new JButton("Select BAM(s)");
        bamButton.addActionListener(this);
        bamPanel.add(bamButton);
        
        bamLocalButton = new JRadioButton("Local");        
        bamLocalButton.setActionCommand("bamLocal");
        bamLocalButton.addActionListener(this);
        bamLocalButton.setSelected(true);
        
        bamRemoteButton = new JRadioButton("Remote");
        bamRemoteButton.setActionCommand("bamRemote");
        bamRemoteButton.addActionListener(this);
        
        ButtonGroup bamFileSource = new ButtonGroup();
        bamFileSource.add(bamLocalButton);
        bamFileSource.add(bamRemoteButton);
        
        bamPanel.add(bamLocalButton);
        bamPanel.add(bamRemoteButton);
        
        contents.add(bamPanel);
        
        JPanel bamListPanel = new JPanel();
        multipleButton = new JRadioButton("");
        multipleButton.setActionCommand("multipleBam");
        multipleButton.addActionListener(this);
        bamListPanel.add(multipleButton);
        bamListPanel.add(new JLabel("BAM list: "));
        bamListField = new JTextField(20);
        bamListField.setEditable(false);
        bamListPanel.add(bamListField);
        bamListButton = new JButton("Select BAM list");
        bamListButton.addActionListener(this);
        bamListButton.setEnabled(false);
        bamListPanel.add(bamListButton);
        
        bamListLocalButton = new JRadioButton("Local");        
        bamListLocalButton.setActionCommand("bamListLocal");
        bamListLocalButton.addActionListener(this);
        bamListLocalButton.setEnabled(false);
        bamListLocalButton.setSelected(true);
        
        bamListRemoteButton = new JRadioButton("Remote");
        bamListRemoteButton.setActionCommand("bamListRemote");
        bamListRemoteButton.addActionListener(this);
        bamListRemoteButton.setEnabled(false);
        
        ButtonGroup bamListSource = new ButtonGroup();
        bamListSource.add(bamListLocalButton);
        bamListSource.add(bamListRemoteButton);
        
        bamListPanel.add(bamListLocalButton);
        bamListPanel.add(bamListRemoteButton);
        
        ButtonGroup bamFileType = new ButtonGroup();
        bamFileType.add(singleButton);
        bamFileType.add(multipleButton);
        
        contents.add(bamListPanel);
        
        
        JPanel varPanel = new JPanel();
        varPanel.add(new JLabel("Variants file: "));
        varField = new JTextField(20);
        varPanel.add(varField);
        varButton = new JButton("Select Variants");
        varButton.addActionListener(this);
        varPanel.add(varButton);
        
        varLocalButton = new JRadioButton("Local");        
        varLocalButton.setActionCommand("varLocal");
        varLocalButton.addActionListener(this);
        varLocalButton.setSelected(true);
        
        varRemoteButton = new JRadioButton("Remote");
        varRemoteButton.setActionCommand("varRemote");
        varRemoteButton.addActionListener(this);
        varRemoteButton.setEnabled(false);
        
        ButtonGroup varFileSource = new ButtonGroup();
        varFileSource.add(varLocalButton);
        varFileSource.add(varRemoteButton);
        
        varPanel.add(varLocalButton);
        varPanel.add(varRemoteButton);
      
        contents.add(varPanel);
        
        JPanel refPanel = new JPanel();
        refPanel.add(new JLabel("Reference file: "));
        refField = new JTextField(20);
        refPanel.add(refField);
        refButton = new JButton("Select Reference");
        refButton.addActionListener(this);
        refPanel.add(refButton);
        
        refLocalButton = new JRadioButton("Local");        
        refLocalButton.setActionCommand("refLocal");
        refLocalButton.addActionListener(this);
        refLocalButton.setSelected(true);
        
        refRemoteButton = new JRadioButton("Remote");
        refRemoteButton.setActionCommand("refRemote");
        refRemoteButton.addActionListener(this);
        refRemoteButton.setEnabled(false);
        
        ButtonGroup refFileSource = new ButtonGroup();
        refFileSource.add(refLocalButton);
        refFileSource.add(refRemoteButton);
        
        refPanel.add(refLocalButton);
        refPanel.add(refRemoteButton);
        
        contents.add(refPanel);
        
        JPanel hostPanel = new JPanel();
        hostPanel.add(new JLabel("Host: "));
        hostField = new JTextField(20);
        hostField.setEnabled(false);
        hostPanel.add(hostField);
        hostPanel.add(new JLabel("Port: "));
        portField = new JTextField("22", 3);
        portField.setEnabled(false);
        hostPanel.add(portField);
        contents.add(hostPanel);
        
        JPanel userPanel = new JPanel();
        userPanel.add(new JLabel("Username: "));
        userField = new JTextField(System.getProperty("user.name"), 10);
        userField.setEnabled(false);
        userPanel.add(userField);
        JPanel passPanel = new JPanel();
        passPanel.add(new JLabel("Password: "));
        pf = new JPasswordField(8);
        pf.setEnabled(false);
        passPanel.add(pf);
        userPanel.add(passPanel);
        contents.add(userPanel);
        
        JPanel butPan = new JPanel();
        JButton okbut = new JButton("OK");
        getRootPane().setDefaultButton(okbut);
        okbut.addActionListener(this);
        butPan.add(okbut);
        JButton cancelbut = new JButton("Cancel");
        cancelbut.addActionListener(this);
        butPan.add(cancelbut);
        contents.add(butPan);

        this.setContentPane(contents);
    }

    public void actionPerformed(ActionEvent e) {
        if (e.getActionCommand().equals("OK")){
            bamFile = bamField.getText();
            bamList = bamListField.getText();
            varFile = varField.getText();
            refFile = refField.getText();
            host = hostField.getText();
            port = portField.getText();
            username = userField.getText();
            password = pf.getPassword();
            filled = true;
            this.dispose();
            
        }else if (e.getActionCommand().equals("Select BAM(s)")) {
//            JFileChooser jfc = new JFileChooser("user.dir");
        	// for testing only
            JFileChooser jfc = new JFileChooser("/Users/jm20/Documents/snug");
            jfc.setDialogTitle("Select BAM(s)");
            jfc.setApproveButtonText("Select BAM(s)");
            jfc.setApproveButtonToolTipText("Select your indexed BAM file or file containing a list of paths to your BAM files");
            
            if (jfc.showOpenDialog(this) == JFileChooser.APPROVE_OPTION){
            	bamField.setText(jfc.getSelectedFile().getAbsolutePath());
            }
            
        } else if (e.getActionCommand().equals("Select Variants")) {
//        	JFileChooser jfc = new JFileChooser("user.dir");
        	// for testing only
        	JFileChooser jfc = new JFileChooser("/Users/jm20/Documents/snug");
        	jfc.setDialogTitle("Select Variant File");
            jfc.setApproveButtonText("Select Variant File");
            jfc.setApproveButtonToolTipText("Select your list of variant sites to visualise");
        	
            if (jfc.showOpenDialog(this) == JFileChooser.APPROVE_OPTION){
                varField.setText(jfc.getSelectedFile().getAbsolutePath());
            }
            
        } else if (e.getActionCommand().equals("Select Reference")) {
//        	JFileChooser jfc = new JFileChooser("user.dir");    
        	// for testing only
        	JFileChooser jfc = new JFileChooser("/Users/jm20/Documents/snug");
        	jfc.setDialogTitle("Select Reference File");
            jfc.setApproveButtonText("Select Reference File");
            jfc.setApproveButtonToolTipText("Select your indexed reference file");

            if (jfc.showOpenDialog(this) == JFileChooser.APPROVE_OPTION){
                refField.setText(jfc.getSelectedFile().getAbsolutePath());
            }
            
        } else if (e.getActionCommand().equals("bamLocal")) {
        	// activate the file browse button
        	bamButton.setEnabled(true);
        
        	if (refLocalButton.isSelected() && varLocalButton.isSelected()) {
        		hostField.setEnabled(false);
            	portField.setEnabled(false);
            	userField.setEnabled(false);
            	pf.setEnabled(false);        		
        	}
        	
        } else if (e.getActionCommand().equals("bamRemote")) {
        	bamButton.setEnabled(false);
        	hostField.setEnabled(true);
        	portField.setEnabled(true);
        	userField.setEnabled(true);
        	pf.setEnabled(true);        
        	
        } else if (e.getActionCommand().equals("bamListLocal")) {
        	// activate the file browse button
        	bamListButton.setEnabled(true);
        
        	if (refLocalButton.isSelected() && varLocalButton.isSelected()) {
        		hostField.setEnabled(false);
            	portField.setEnabled(false);
            	userField.setEnabled(false);
            	pf.setEnabled(false);        		
        	}
        	
        } else if (e.getActionCommand().equals("bamListRemote")) {
        	bamListButton.setEnabled(false);
        	hostField.setEnabled(true);
        	portField.setEnabled(true);
        	userField.setEnabled(true);
        	pf.setEnabled(true);        
        	
        } else if (e.getActionCommand().equals("varLocal")) {
        	varButton.setEnabled(true);        	
        	if (refLocalButton.isSelected() && bamLocalButton.isSelected()) {
        		hostField.setEnabled(false);
            	portField.setEnabled(false);
            	userField.setEnabled(false);
            	pf.setEnabled(false);        		
        	}
        	
        } else if (e.getActionCommand().equals("varRemote")) {
        	varButton.setEnabled(false);
        	hostField.setEnabled(true);
        	portField.setEnabled(true);
        	userField.setEnabled(true);
        	pf.setEnabled(true);
        	
        } else if (e.getActionCommand().equals("refLocal")) {
        	refButton.setEnabled(true);
        	if (varLocalButton.isSelected() && bamLocalButton.isSelected()) {
        		hostField.setEnabled(false);
            	portField.setEnabled(false);
            	userField.setEnabled(false);
            	pf.setEnabled(false);        		
        	}
        	
        } else if (e.getActionCommand().equals("refRemote")) {
        	refButton.setEnabled(false);
        	hostField.setEnabled(true);
        	portField.setEnabled(true);
        	userField.setEnabled(true);
        	pf.setEnabled(true);
        	
        } else if (e.getActionCommand().equals("singleBam")) {
            bamField.setEditable(true);
            bamButton.setEnabled(true);
            bamLocalButton.setEnabled(true);
            bamRemoteButton.setEnabled(true);
            bamListField.setEditable(false);
            bamListButton.setEnabled(false);
            bamListLocalButton.setEnabled(false);
            bamListRemoteButton.setEnabled(false);        	
        } else if (e.getActionCommand().equals("multipleBam")) {
        	bamField.setEditable(false);
            bamButton.setEnabled(false);
            bamLocalButton.setEnabled(false);
            bamRemoteButton.setEnabled(false);
            bamListField.setEditable(true);
            bamListButton.setEnabled(true);
            bamListLocalButton.setEnabled(true);
            bamListRemoteButton.setEnabled(true);
        } else if (e.getActionCommand().equals("Cancel")){
        	this.dispose();
        }
    }

	public String getRef() {
		return refFile;
	}

	public boolean localRef() {
		return refLocalButton.isSelected();
	}
	
	public boolean remoteRef() {
		return refRemoteButton.isSelected();
	}
	
	public String getVar() {
		return varFile;
	}

	public boolean localVar() {
		return varLocalButton.isSelected();
	}
	
	public boolean remoteVar() {
		return varRemoteButton.isSelected();
	}
	
	public String getBam() {
		return bamFile;
	}
	
	public String getBamList() {
		return bamList;
	}
	
	public boolean localBam() {
		return bamLocalButton.isSelected();
	}
	
	public boolean remoteBam() {
		return bamRemoteButton.isSelected();
	}
	
	public String getHost() {
		return host;
	}

	public int getPort() {
		return Integer.parseInt(port);
	}

	public String getUsername() {
		return username;
	}

	public String getPassword() {
		return new String(password);
	}
	
	public boolean isFilled() {
		return filled;
	}
}