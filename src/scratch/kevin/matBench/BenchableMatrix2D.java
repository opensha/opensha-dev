package scratch.kevin.matBench;

public interface BenchableMatrix2D {
	
	public int rows();
	public int cols();
	public double get(int row, int col);
	public void set(int row, int col, double value);
	
	public BenchableMatrix2D mult(BenchableMatrix2D mat);
	
	public long getMultTime();
	public void resetMultTime();

}
