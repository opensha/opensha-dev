package scratch.kevin;

import java.util.ArrayDeque;
import java.util.Collections;
import java.util.Deque;
import java.util.List;

import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.primitives.Ints;

public class PrisonerTest {

	public static void main(String[] args) {
		int[] positions = new int[100];
		for (int i=0; i<positions.length; i++)
			positions[i] = i;
		
		int numTrials = 1000000;
		
		int numSucceeded = 0;
		
		for (int t=0; t<numTrials; t++) {
			List<Integer> myPos = Ints.asList(positions);
			Collections.shuffle(myPos);
			
			boolean print = t % 100000 == 0;
			
			boolean trialSuccess = true;
//			for (int p=0; p<2; p++) {
////			for (int p=0; p<positions.length; p++) {
//				int startIndex, endIndex;
////				if (p % 2 == 0) {
//					startIndex = 0;
//					endIndex = positions.length/2;
////				} else {
////					startIndex = positions.length/2;
////					endIndex = positions.length;
////				}
////				int startIndex = p;
////				int endIndex = p+positions.length/2;
////				int startIndex = 0;
////				int endIndex = positions.length;
//				
//				boolean success = false;
//				List<String> myTries = Lists.newArrayList();
//				int myTry = 0;
//				boolean up = true;
//				if (p % 2 == 1) {
//					myTry = 99;
//					up = false;
//				}
////				if (p == 2)
////					myTry = 50;
//				for (int i=startIndex; i<endIndex; i++) {
////					int index = i % positions.length;
////					myTry = (p + 2*i) % positions.length;
//					int val = myPos.get(myTry);
//					if (val == p) {
//						success = true;
//						if (!print)
//							break;
//					}
//					if (val == p)
//						myTries.add("*"+myTry+"*");
//					else
//						myTries.add(myTry+"");
////					if (p % 2 == 0) {
////						if (val == p+1)
////							myTry = 50;
////						else
////							myTry++;
////					} else {
////						if (val == p-1)
////							myTry = 0;
////						else if (myTry >= 50)
////							myTry--;
////						else
////							myTry++;
////					}
//					
////					myTry++;
//					if (p % 2 == 0) {
//						if (val == p+1) {
//							myTry = 100;
//							up = false;
//						}
//					} else {
//						if (val == p-1) {
//							myTry = -1;
//							up = true;
//						}
//					}
//					if (up)
//						myTry++;
//					else
//						myTry--;
//				}
//				if (!success)
//					Preconditions.checkState(myTries.size() == 50);
//				if (print) {
//					if (p == 0)
//						System.out.println("Trial "+t);
//					System.out.println("\t"+p+": "+Joiner.on(",").join(myTries));
//				}
//				if (!success) {
//					trialSuccess = false;
//					if (!print)
//						break;
//				}
//			}
			
//			for (int p=0; p<positions.length; p++) {
			for (int p=0; p<positions.length; p++) {
				boolean success = false;
				
				List<Integer> trials = Lists.newArrayList();
				
//				int startIndex = p % 2;
//				int[] startingPoints = { 0, 1, 0, 98, 99, 0 };
//				boolean[] ups = { true, true, true, false, false, true };
				
//				int curTest = startingPoints[startIndex];
//				boolean up = ups[startIndex];
				
				int curTest = p;
				
				for (int trial=0; trial<50; trial++) {
					int val = myPos.get(curTest);
					
					trials.add(curTest);
					
					if (val == p) {
						success = true;
						if (!print)
							break;
					}
					
					curTest = val;
				}
				
//				Deque<Integer> evensLeft = new ArrayDeque<Integer>();
//				Deque<Integer> oddsLeft = new ArrayDeque<Integer>();
//				for (int i=0; i<positions.length; i++) {
//					if (i % 2 == 0)
//						evensLeft.add(i);
//					else
//						oddsLeft.add(i);
//				}
//				
//				boolean even = p % 2 == 0;
//				
//				for (int trial=0; trial<50; trial++) {
//					int curTest;
//					if (even)
//						curTest = evensLeft.pop();
//					else
//						curTest = oddsLeft.pop();
//					int val = myPos.get(curTest);
//					
//					trials.add(curTest);
//					
//					if (val == p) {
//						success = true;
////						break;
//					}
////					if (val == p + 1 || val == p - 1) {
////					if (val == p - 1 || val == p + 1) {
////						startIndex++;
////						curTest = startingPoints[startIndex];
////						up = ups[startIndex];
////					} else {
////						if (up)
////							curTest += 2;
////						else
////							curTest -= 2;
////					}
////					if (val <= p+1)
////						curTest += 1;
////					else
////						curTest += 2;
//					
//					if (val == p+1 || val == p -1)
//						even = !even;
//				}
				
				if (print) {
					if (p == 0)
						System.out.println("Trial "+t);
					System.out.println("\t"+p+": "+Joiner.on(",").join(trials));
				}
				
				if (!success) {
					trialSuccess = false;
//					break;
				}
			}
			
			if (trialSuccess)
				numSucceeded++;
		}
		
		double percent = 100d*numSucceeded/(double)numTrials;
		System.out.println(numSucceeded+"/"+numTrials+" succeeded ("+(float)percent+" %)");
	}

}
