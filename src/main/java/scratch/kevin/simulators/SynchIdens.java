package scratch.kevin.simulators;

import java.awt.Color;
import java.util.List;

import org.opensha.sha.simulators.iden.ElementMagRangeDescription;
import org.opensha.sha.simulators.iden.LogicalOrRupIden;
import org.opensha.sha.simulators.iden.RuptureIdentifier;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

public class SynchIdens {
	
	public static RuptureIdentifier build(double minMag, double maxMag, SynchFaults... faults) {
		return build(null, minMag, maxMag, faults);
	}
	
	public static RuptureIdentifier build(String name, double minMag, double maxMag, SynchFaults... faults) {
		String magStr;
		if (maxMag >= 10)
			magStr = SimAnalysisCatLoader.getMagStr(minMag)+"+";
		else
			magStr = SimAnalysisCatLoader.getMagStr(minMag)+"=>"+SimAnalysisCatLoader.getMagStr(maxMag);
		Preconditions.checkArgument(faults.length > 0, "Must supply at least one fault");
		if (faults.length == 1)
			if (name == null)
				return new ElementMagRangeDescription(faults[0].name+" "+magStr, faults[0].allcal2ID, minMag, maxMag);
			else
				return new ElementMagRangeDescription(name+" "+magStr, faults[0].allcal2ID, minMag, maxMag);
		
		// multiple - do logical or
		List<RuptureIdentifier> idens = Lists.newArrayList();
		for (SynchFaults fault : faults)
			idens.add(new ElementMagRangeDescription(fault.name+" "+magStr, fault.allcal2ID, minMag, maxMag));
		LogicalOrRupIden iden = new LogicalOrRupIden(idens);
		if (name != null)
			iden.setName(name+" "+magStr);
		return iden;
	}
	
	public static List<RuptureIdentifier> getIndividualFaults(double minMag, double maxMag, SynchFaults... faults) {
		List<RuptureIdentifier> idens = Lists.newArrayList();
		for (SynchFaults fault : faults)
			idens.add(build(minMag, maxMag, fault));
		return idens;
	}
	
	public static List<Color> getStandardColors() {
		return Lists.newArrayList(Color.BLACK, Color.BLUE, Color.RED, Color.MAGENTA,
				Color.GREEN, Color.CYAN, Color.GRAY, Color.ORANGE, Color.YELLOW, Color.PINK);
	}
	
	public static List<RuptureIdentifier> getStandardSoCal() {
		List<RuptureIdentifier> idens = getIndividualFaults(7d, 10d,
				SynchFaults.SAF_MOJAVE, SynchFaults.SAF_CARRIZO,
				SynchFaults.SAF_COACHELLA, SynchFaults.SAF_CHOLAME,
				SynchFaults.GARLOCK_WEST, SynchFaults.SAN_JACINTO);
		
		// TODO now add multi basin
		
		return idens;
	}
	
	public static List<RuptureIdentifier> getStandardNorCal() {
		return getIndividualFaults(7d, 10d,
				SynchFaults.SAF_SANTA_CURZ, SynchFaults.SAF_MID_PENINSULA,
				SynchFaults.CALAVERAS, SynchFaults.HAYWARD,
				SynchFaults.SAN_GREGORIO, SynchFaults.RODGERS_CREEK,
				SynchFaults.MAACAMA_CENTRAL, SynchFaults.SAF_N_COAST_ON);
	}
	
	public static enum SynchFaults {
		//						NAME							ALLCAL2 ID
		SAF_MOJAVE(				"SAF Mojave",					1246),
		SAF_CARRIZO(			"SAF Carrizo",					1026),
		SAF_COACHELLA(			"SAF Coachella",				1602),
		SAF_CHOLAME(			"SAF Cholame",					966),
		GARLOCK_WEST(			"Garlock",						6036),
		SAN_JACINTO(			"San Jacinto",					1931),
		ELSINORE(				"Elsinore",						2460),
		PUENTE_HILLS(			"Puente Hills",					11829),
		NEWP_INGL_ONSHORE(		"Newport-Inglewood (onshore)",	7672),
		SAF_SANTA_CURZ(			"SAF Santa Cruz",				663),
		SAF_MID_PENINSULA(		"SAF Mid Peninsula",			529),
		CALAVERAS(				"Calaveras",					3692),
		HAYWARD(				"Hayward",						3334),
		SAN_GREGORIO(			"San Gregorio",					2782),
		RODGERS_CREEK(			"Rodgers Creek",				3238),
		MAACAMA_CENTRAL(		"Maacama Central",				3098),
		SAF_N_COAST_ON(			"SAF N Coast (onshore)",		389);
		
		private String name;
		private int allcal2ID;
		
		private SynchFaults(String name, int allcal2ID) {
			this.name = name;
			this.allcal2ID = allcal2ID;
		}
	}

}
