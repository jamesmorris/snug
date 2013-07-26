
public class Variant {

	String chr;
	int pos;
	String name;
	
	public Variant(String c, int p) {
		this.chr = c;
		this.pos = p;
		this.name = c + ":" + Integer.toString(p); 
	}

	public String getChr() {
		return chr;
	}

	public void setChr(String chr) {
		this.chr = chr;
	}

	public int getPos() {
		return pos;
	}

	public void setPos(int pos) {
		this.pos = pos;
	}

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

}
