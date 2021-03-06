import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;



public class VariantListFile {

	ArrayList<Variant> variants = new ArrayList<Variant>();
	
	public VariantListFile(String filename) throws IOException, UnrecognisedVariantFileFormat {
		
		BufferedReader listReader = new BufferedReader(new FileReader(filename));
		
		String currentline;
		
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
					String chr = bits[0];
					int pos = PosStr2Int(bits[1]);
					variants.add(new Variant(chr, pos));	
				} else {
					variants = null;
					throw new  UnrecognisedVariantFileFormat("Can't parse variant chromosome and position '" + currentline + "'");
					
				}					
			}
						
		}
		listReader.close();
	}

	public ArrayList<Variant> getVariants() {
		return variants;
	}

	public void setVariants(ArrayList<Variant> v) {
		this.variants = v;
	}
	
    private int PosStr2Int(String string) {        
        return Integer.parseInt(string);
    }
	
}
