import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.Iterator;

import org.junit.Test;


public class snugViewerTest {

	@Test
	public void testOpenVarientFile() {		
		VariantListFile variantList = new VariantListFile("/Users/jm20/Documents/snug/test_variant_file");
		
		ArrayList<Variant> variants = variantList.getVariants();
		
		assertEquals("number of variants in test_variant_file", 11, variants.size());
		
		Iterator<Variant> variantItr = variants.iterator();
		
		Variant v = variantItr.next();		
		assertEquals("check parse ':' delimeter chromosome value", "20", v.getChr());
		assertEquals("check parse ':' delimeter position value", 14370, v.getPos());
		
		v = variantItr.next();
		assertEquals("check parse X chromosome", "23", v.getChr());
		
		v = variantItr.next();
		assertEquals("check parse Y chromosome", "24", v.getChr());
				
		v = variantItr.next();
		assertEquals("check parse XY chromosome", "26", v.getChr());
				
		v = variantItr.next();
		assertEquals("check parse MT chromosome", "25", v.getChr());
				
		v = variantItr.next();
		assertEquals("check parse chr prefix", "20", v.getChr());
				
		v = variantItr.next();
		assertEquals("check parse Chr prefix", "20", v.getChr());
		
		v = variantItr.next();
		assertEquals("check parse CHR prefix", "20", v.getChr());
	
		v = variantItr.next();		
		assertEquals("check parse '\t' delimeter chromosome value", 20, v.getChr());
		assertEquals("check parse '\t' delimeter chromosome value", 14370, v.getPos());
		
		v = variantItr.next();		
		assertEquals("check parse ' ' delimeter chromosome value", 20, v.getChr());
		assertEquals("check parse ' ' delimeter chromosome value", 14370, v.getPos());			
		
		v = variantItr.next();		
		assertEquals("check parse ',' delimeter chromosome value", 20, v.getChr());
		assertEquals("check parse ',' delimeter chromosome value", 14370, v.getPos());
		
	}

}
