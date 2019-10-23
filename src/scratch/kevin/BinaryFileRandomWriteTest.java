package scratch.kevin;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FilterOutputStream;
import java.io.IOException;
import java.io.RandomAccessFile;

public class BinaryFileRandomWriteTest {
	
	private static void printFile(File binaryFile) throws IOException {
		DataInputStream in = new DataInputStream(new FileInputStream(binaryFile));
		
		System.out.println("CONTENTS");
		int index = 0;
		while (true) {
			try {
				System.out.println("\t"+(index++)+": "+in.readInt());
			} catch (IOException e) {
				System.out.println("---END---");
				break;
			}
		}
		in.close();
	}

	public static void main(String[] args) throws IOException {
		File file = new File("/tmp/bin_test.bin");
		
		DataOutputStream out = new DataOutputStream(new FileOutputStream(file));
		for (int i=0; i<10; i++)
			out.writeInt(i);
		
		out.close();
		
		printFile(file);
		
		FileOutputStream appendOut = new FileOutputStream(file, true);
		
		out = new DataOutputStream(appendOut);
		for (int i=0; i<10; i++)
			out.writeInt(i);
		
		out.close();
		
		printFile(file);
		
		// now try to rewrite the last part
		
		RandomAccessFile raFile = new RandomAccessFile(file, "rw");
		raFile.seek(40l);
		
		appendOut = new FileOutputStream(raFile.getFD());
		
		out = new DataOutputStream(appendOut);
		for (int i=0; i<20; i++)
			out.writeInt(i+1000);
		
		out.close();
		
		printFile(file);
	}

}
