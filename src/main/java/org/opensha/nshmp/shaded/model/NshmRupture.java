package org.opensha.nshmp.shaded.model;

import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.faultSurface.RuptureSurface;

import org.opensha.nshmp.shaded.NshmpMaths;

class NshmRupture extends ProbEqkRupture {

	final double rate;
	final double weight;

	NshmRupture(
			double mag,
			double rake,
			double rate,
			double weight,
			double duration,
			RuptureSurface surface) {

		super(
				mag, rake,
				NshmpMaths.rateToProbability(rate * weight, duration),
				surface, null);

		this.rate = rate;
		this.weight = weight;
	}
}
