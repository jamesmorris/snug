
public class BasePosition {
	int a = 0;
	int t = 0;
	int g = 0;
	int c = 0;
	
	public int getA() {
		return a;
	}

	public void setA(int a) {
		this.a = a;
	}

	public int getT() {
		return t;
	}

	public void setT(int t) {
		this.t = t;
	}

	public int getG() {
		return g;
	}

	public void setG(int g) {
		this.g = g;
	}

	public int getC() {
		return c;
	}

	public void setC(int c) {
		this.c = c;
	}

	public BasePosition() {
		
	}

	public void addBase(String base) {
		if (base.equals("A")) {
			setA(getA() + 1);
		} else if (base.equals("T")) {
			setT(getT() + 1);
		} else if (base.equals("C")) {
			setC(getC() + 1);
		} else if (base.equals("G")) {
			setG(getG() + 1);
		} 
	}
}
