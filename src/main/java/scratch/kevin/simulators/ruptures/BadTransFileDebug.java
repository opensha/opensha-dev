package scratch.kevin.simulators.ruptures;

import java.io.File;
import java.io.IOException;
import java.nio.ByteOrder;

import org.opensha.sha.simulators.srf.RSQSimStateTime;
import org.opensha.sha.simulators.srf.RSQSimStateTransitionFileReader;
import org.opensha.sha.simulators.srf.RSQSimStateTransitionFileReader.TransVersion;

import com.google.common.base.Preconditions;

public class BadTransFileDebug {

	public static void main(String[] args) throws IOException {
		File dir = new File("/project/scec_608/rsqsim/catalogs/shaw/rundir5672/");
		File transFile = new File(dir, "trans..out");
		
		TransVersion version = TransVersion.CONSOLIDATED_RELATIVE;
		ByteOrder order = ByteOrder.LITTLE_ENDIAN;
		
		RSQSimStateTransitionFileReader transRead = new RSQSimStateTransitionFileReader(
				transFile, order, version);
		
		long numVals = transFile.length() / version.bytesPerRecord;
		System.out.println("Have "+numVals+" possible transitions");
		
//		long startIndex = 925244444l;
		long startIndex = 900000000l;
		int bufferLen = 100000;
		
		long firstBadIndex = -1;
		long lastBadIndex = -1;
		long numBad = 0;
		long consecutiveBad = 0;
		long consecutiveGood = 0;
		RSQSimStateTime lastGoodTrans = null;
		RSQSimStateTime firstBadTrans = null;
		RSQSimStateTime lastBadTrans = null;
		int curEventID = -1;
		long firstIndexForEvent = -1;
		
		System.out.println("Starting read at index "+startIndex);
		
		long index = startIndex;
		while (index < numVals) {
			RSQSimStateTime[] trans = transRead.read(index, bufferLen);
			for (int i=0; i<trans.length; i++) {
				if (trans[i].eventID > 0 && trans[i].eventID != curEventID) {
					firstIndexForEvent = index;
					curEventID = trans[i].eventID;
				}
				if (trans[i].patchID <= 0) {
					if (firstBadIndex < 0) {
						System.out.println("First bad transition found at index "+index+":\t"+trans[i]);
						System.out.println("Prior good transition:\t"+lastGoodTrans);
						System.out.println("First index for this event ("+curEventID+"):\t"
								+firstIndexForEvent+" at position "+(firstIndexForEvent*version.bytesPerRecord));
						firstBadTrans = trans[i];
						firstBadIndex = index;
					} else if (consecutiveGood > 0l) {
						System.out.println("We're bad again at "+index+" after "+consecutiveGood+" consecutive good:\t"+trans[i]);
						consecutiveBad = 0l;
					}
					consecutiveGood = 0l;
					consecutiveBad++;
					numBad++;
					lastBadIndex = index;
					lastBadTrans = trans[i];
				} else {
					if (consecutiveBad > 0l) {
						System.out.println("We're good again at "+index+" after "+consecutiveBad+" consecutive bad:\t"+trans[i]);
						consecutiveGood = 0l;
					}
					consecutiveBad = 0l;
					consecutiveGood++;
					lastGoodTrans = trans[i];
				}
				index++;
			}
		}
		
		System.out.println("Done, ended with a streak of "+consecutiveBad+" bad/"+consecutiveGood+" good");
		System.out.println("Last good before bad at index "+(firstBadIndex-1)+":\t"+lastGoodTrans);
		System.out.println("First bad at index "+firstBadIndex+":\t"+firstBadTrans);
		System.out.println("Last bad at index "+lastBadIndex+":\t"+lastBadTrans);
		System.out.println("Detected "+numBad+" bad transitions");
	}

}
