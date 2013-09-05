import javax.swing.*;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;

public class DataConnectionDialog extends JDialog implements ActionListener {
    private String bamFile;
    private String varFile;
    private String refFile;
    private JTextField bamField;
    private JTextField varField;
    private JTextField refField;
    private JCheckBox emptyIt;

    public DataConnectionDialog(JFrame parent){
        super(parent,"Data Connection",true);

        JPanel contents = new JPanel();
        contents.setLayout(new BoxLayout(contents,BoxLayout.Y_AXIS));

        JPanel bamPanel = new JPanel();
        bamPanel.add(new JLabel("BAM file: "));
        bamField = new JTextField(20);
        bamPanel.add(bamField);
        JButton bamButton = new JButton("Select BAM(s)");
        bamButton.addActionListener(this);
        bamPanel.add(bamButton);
        contents.add(bamPanel);
        
        JPanel varPanel = new JPanel();
        varPanel.add(new JLabel("Variants file: "));
        varField = new JTextField(20);
        varPanel.add(varField);
        JButton varButton = new JButton("Select Variants");
        varButton.addActionListener(this);
        varPanel.add(varButton);
        contents.add(varPanel);
        
        JPanel refPanel = new JPanel();
        refPanel.add(new JLabel("Reference file: "));
        refField = new JTextField(20);
        refPanel.add(refField);
        JButton refButton = new JButton("Select Reference");
        refButton.addActionListener(this);
        refPanel.add(refButton);
        contents.add(refPanel);
        
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
            varFile = varField.getText();
            refFile = refField.getText();
            
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
        } else if (e.getActionCommand().equals("Cancel")){
            this.dispose();
        }
    }

	public String getRef() {
		return refFile;
	}

	public String getVar() {
		return varFile;
	}

	public String getBam() {
		return bamFile;
	}
}