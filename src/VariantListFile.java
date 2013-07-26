import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedList;


public class VariantListFile {

	ArrayList<Variant> variants = new ArrayList<Variant>();
	
	public VariantListFile(String filename) {
		BufferedReader listReader = null;
		try {
			listReader = new BufferedReader(new FileReader(filename));
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
					String[] bits;					
					if (currentline.split(":").length == 2) {
						bits = currentline.split(":");
					} else if (currentline.split("\t").length == 2) {
						bits = currentline.split("\t");
					} else if (currentline.split(" +").length == 2) {
						bits = currentline.split(" +");
					} else if  (currentline.split(",").length == 2) {
						bits = currentline.split(",");
					} else if (currentline.split("-").length == 2) {
						bits = currentline.split("-");
					} else {
						bits = null;
					}					
					if (bits != null) {						
						int[] variant = new int[2];
						String chr = bits[0];
						int pos = PosStr2Int(bits[1]);
						variants.add(new Variant(chr, pos));	
					} else {
						variants = null;
						System.out.println("Can't parse variant chromosome and position '" + currentline + "'");
					}					
				}
						
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}

	public ArrayList<Variant> getVariants() {
		return variants;
	}

	public void setVariants(ArrayList<Variant> v) {
		this.variants = v;
	}
	
    private int PosStr2Int(String string) {
        int pos = 0;
        try {
            pos = Integer.parseInt(string);
        } catch (NumberFormatException nfe) {
            System.out.println("position is not a number '" + string + "'");
        }
        return pos;
    }
	
}
