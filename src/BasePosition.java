
public class BasePosition {
	int ap = 0;
	int tp = 0;
	int gp = 0;
	int cp = 0;
	int am = 0;
	int tm = 0;
	int gm = 0;
	int cm = 0;
	
	public int getAplus() {
		return ap;
	}

	public void setAplus(int a) {
		this.ap = a;
	}

	public int getTplus() {
		return tp;
	}

	public void setTplus(int t) {
		this.tp = t;
	}

	public int getGplus() {
		return gp;
	}

	public void setGplus(int g) {
		this.gp = g;
	}

	public int getCplus() {
		return cp;
	}

	public void setCplus(int c) {
		this.cp = c;
	}

	public BasePosition() {
		
	}
	public int getAminus() {
		return am;
	}

	public void setAminus(int a) {
		this.am = a;
	}

	public int getTminus() {
		return tm;
	}

	public void setTminus(int t) {
		this.tm = t;
	}

	public int getGminus() {
		return gm;
	}

	public void setGminus(int g) {
		this.gm = g;
	}

	public int getCminus() {
		return cm;
	}

	public void setCminus(int c) {
		this.cm = c;
	}
	
	public void addBase(String base) {
		if (base.equals("+A")) {
			setAplus(getAplus() + 1);
		} else if (base.equals("+T")) {
			setTplus(getTplus() + 1);
		} else if (base.equals("+C")) {
			setCplus(getCplus() + 1);
		} else if (base.equals("+G")) {
			setGplus(getGplus() + 1);
		} else if (base.equals("-A")) {
			setAminus(getAminus() + 1);
		} else if (base.equals("-T")) {
			setTminus(getTminus() + 1);
		} else if (base.equals("-C")) {
			setCminus(getCminus() + 1);
		} else if (base.equals("-G")) {
			setGminus(getGminus() + 1);
		}  
	}
}
