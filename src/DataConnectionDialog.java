import javax.swing.*;

import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;

public class DataConnectionDialog extends JDialog implements ActionListener {
    private String bamFile;
    private String bamList;
    private String bamVCFPath;
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
	JRadioButton localFilesButton;
	JRadioButton remoteFilesButton;
	JRadioButton singleButton;
    JRadioButton multipleButton;
	private boolean filled;
        
    JButton bamButton;
    JButton bamListButton;
    JButton varButton;
    JButton refButton;
	private JRadioButton vcfButton;
	private JTextField vcfbamField;
	
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
        bamPanel.add(new JLabel("Single BAM: "));
        bamField = new JTextField(20);
        bamPanel.add(bamField);
        bamButton = new JButton("Select BAM");
        bamButton.addActionListener(this);
        bamPanel.add(bamButton);
        
        contents.add(bamPanel);
        
        JPanel bamListPanel = new JPanel();
        multipleButton = new JRadioButton("");
        multipleButton.setActionCommand("multipleBam");
        multipleButton.addActionListener(this);
        bamListPanel.add(multipleButton);
        bamListPanel.add(new JLabel("Multiple BAMs: "));
        bamListField = new JTextField(20);
        bamListField.setEditable(false);
        bamListPanel.add(bamListField);
        bamListButton = new JButton("Select BAM list");
        bamListButton.addActionListener(this);
        bamListButton.setEnabled(false);
        bamListPanel.add(bamListButton);
        
        contents.add(bamListPanel);
        
        JPanel bamSourcePanel = new JPanel();
        
        localFilesButton = new JRadioButton("Local file(s)");        
        localFilesButton.setActionCommand("localFiles");
        localFilesButton.addActionListener(this);
        localFilesButton.setSelected(true);
        
        remoteFilesButton = new JRadioButton("Remote file(s)");
        remoteFilesButton.setActionCommand("remoteFiles");
        remoteFilesButton.addActionListener(this);
        
        ButtonGroup bamListSource = new ButtonGroup();
        bamListSource.add(localFilesButton);
        bamListSource.add(remoteFilesButton);
        
        bamSourcePanel.add(localFilesButton);
        bamSourcePanel.add(remoteFilesButton);
        
        contents.add(bamSourcePanel);
        
//        JPanel bamVCFPanel = new JPanel();
//        vcfButton = new JRadioButton("");
//        vcfButton.setActionCommand("vcfBam");
//        vcfButton.addActionListener(this);
//        bamVCFPanel.add(vcfButton);
//        bamVCFPanel.add(new JLabel("Load BAMs From VCF: "));
//        vcfbamField = new JTextField(20);
//        vcfbamField.setEditable(false);
//        bamVCFPanel.add(vcfbamField);
//        
//        contents.add(bamVCFPanel);
        
        // assume that the bams directory is remote this is quite a safe assumption because a directory containing 100s of sample bams will be very large
        
        ButtonGroup bamFileType = new ButtonGroup();
        bamFileType.add(singleButton);
        bamFileType.add(multipleButton);
//        bamFileType.add(vcfButton);
        
        JPanel varPanel = new JPanel();
        varPanel.add(new JLabel("Variants file: "));
        varField = new JTextField(20);
        varPanel.add(varField);
        varButton = new JButton("Select Variants");
        varButton.addActionListener(this);
        varPanel.add(varButton);
 
        contents.add(varPanel);
        
        JPanel refPanel = new JPanel();
        refPanel.add(new JLabel("Reference file: "));
        refField = new JTextField(20);
        refPanel.add(refField);
        refButton = new JButton("Select Reference");
        refButton.addActionListener(this);
        refPanel.add(refButton);
        
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
//            bamVCFPath = vcfbamField.getText();
            varFile = varField.getText();
            refFile = refField.getText();
            host = hostField.getText();
            port = portField.getText();
            username = userField.getText();
            password = pf.getPassword();
            filled = true;
            this.dispose();
            
        } else if (e.getActionCommand().equals("Select BAM")) {
            JFileChooser jfc = new JFileChooser("user.dir");
            jfc.setDialogTitle("Select BAM(s)");
            jfc.setApproveButtonText("Select BAM(s)");
            jfc.setApproveButtonToolTipText("Select your indexed BAM file");
            
            if (jfc.showOpenDialog(this) == JFileChooser.APPROVE_OPTION){
            	bamField.setText(jfc.getSelectedFile().getAbsolutePath());
            }
            
        } else if (e.getActionCommand().equals("Select BAM list")) {
          JFileChooser jfc = new JFileChooser("user.dir");
          jfc.setDialogTitle("Select BAM list");
          jfc.setApproveButtonText("Select BAM list");
          jfc.setApproveButtonToolTipText("Select a file containing a list of paths to your BAM files");
          
          if (jfc.showOpenDialog(this) == JFileChooser.APPROVE_OPTION){
          	bamListField.setText(jfc.getSelectedFile().getAbsolutePath());
          }
          
        } else if (e.getActionCommand().equals("Select Variants")) {
        	JFileChooser jfc = new JFileChooser("user.dir");
        	jfc.setDialogTitle("Select Variant File");
            jfc.setApproveButtonText("Select Variant File");
            jfc.setApproveButtonToolTipText("Select your list of variant sites to visualise");
        	
            if (jfc.showOpenDialog(this) == JFileChooser.APPROVE_OPTION){
                varField.setText(jfc.getSelectedFile().getAbsolutePath());
            }
            
        } else if (e.getActionCommand().equals("Select Reference")) {
        	JFileChooser jfc = new JFileChooser("user.dir");    
        	jfc.setDialogTitle("Select Reference File");
            jfc.setApproveButtonText("Select Reference File");
            jfc.setApproveButtonToolTipText("Select your indexed reference file");

            if (jfc.showOpenDialog(this) == JFileChooser.APPROVE_OPTION){
                refField.setText(jfc.getSelectedFile().getAbsolutePath());
            }
            
        } else if (e.getActionCommand().equals("localFiles")) {        
        	hostField.setEnabled(false);
            portField.setEnabled(false);
            userField.setEnabled(false);
            pf.setEnabled(false);        		
        } else if (e.getActionCommand().equals("remoteFiles")) {
        	hostField.setEnabled(true);
        	portField.setEnabled(true);
        	userField.setEnabled(true);
        	pf.setEnabled(true);        
        } else if (e.getActionCommand().equals("singleBam")) {
            bamField.setEditable(true);
            bamButton.setEnabled(true);
            bamListField.setEditable(false);
            bamListButton.setEnabled(false);
//            vcfbamField.setEditable(false);
        } else if (e.getActionCommand().equals("multipleBam")) {
        	bamField.setEditable(false);
            bamButton.setEnabled(false);
            bamListField.setEditable(true);
            bamListButton.setEnabled(true);
//            vcfbamField.setEditable(false);
        } else if (e.getActionCommand().equals("vcfBam")) {
        	bamField.setEditable(false);
            bamButton.setEnabled(false);
            bamListField.setEditable(false);
            bamListButton.setEnabled(false);            
//            vcfbamField.setEditable(true);            
        } else if (e.getActionCommand().equals("Cancel")){
        	this.dispose();
        }
    }

	public String getRef() {
		return refFile.trim();
	}
	
	public String getVar() {
		return varFile.trim();
	}
	
	public String getBam() {
		return bamFile.trim();
	}
	
	public String getBamList() {
		return bamList.trim();
	}
	
	public String getVCFBamPath() {
		return bamVCFPath;
	}
	
	public boolean localBam() {
		return localFilesButton.isSelected();
	}
	
	public boolean remoteBam() {
		return remoteFilesButton.isSelected();
	}
	
	public boolean singleBam() {
		return singleButton.isSelected();
	}
	
	public boolean multipleBam() {
		return multipleButton.isSelected();
	}
	
	public boolean vcfBam() {
		return vcfButton.isSelected();
	}
	
	public String getHost() {
		return host.trim();
	}

	public int getPort() {
		return Integer.parseInt(port);
	}

	public String getUsername() {
		return username.trim();
	}

	public String getPassword() {
		return new String(password);
	}
	
	public boolean isFilled() {
		return filled;
	}
}