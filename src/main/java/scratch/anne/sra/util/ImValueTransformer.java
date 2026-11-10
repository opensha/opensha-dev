package scratch.anne.sra.util;

/**
 * Strategy interface for transforming intensity measure (IM) arrays and associated
 * response arrays in {@link SimpleImResponseInterface} preparation.
 *
 * <p>This allows flexible interpolation, resampling, or other modifications of the
 * IM values before they are stored in a {@link SimpleImResponseLibrary}.
 */
public interface ImValueTransformer {

    /**
     * Encapsulates the transformed arrays.
     */
    class Result {
        /** Transformed intensity measure values (IM bin edges) */
        public final double[] imValues;
        /** Transformed response values (same length as imValues) */
        public final double[] respValues;

        public Result(double[] imValues, double[] respValues) {
            this.imValues = imValues;
            this.respValues = respValues;
        }
    }

    /**
     * Transform the provided IM and response arrays into new arrays.
     *
     * @param imValues original intensity measure values
     * @param respValues original response values (damage ratios)
     * @return transformed IM and response arrays encapsulated in a Result
     */
    Result transform(double[] imValues, double[] respValues);

    // ------------------- Built-in implementations -------------------

    /**
     * Identity transformer that returns the original arrays unchanged.
     * Use this when no transformation (interpolation, resampling) is needed.
     */
    class NoImTransformation implements ImValueTransformer {
        @Override
        public Result transform(double[] imValues, double[] respValues) {
            // simply return clones to prevent accidental mutation
            return new Result(imValues.clone(), respValues.clone());
        }
    }

    /**
     * Transformer that resamples IM and response arrays onto a log-spaced IM grid.
     * Linear interpolation is used for the response values.
     */
    class LogInterpTransformation implements ImValueTransformer {

        private final double logStep; // spacing in log10(IM)

        /**
         * @param logStep spacing in log10(IM), e.g., 0.025 gives ~90 bins per decade
         */
        public LogInterpTransformation(double logStep) {
            if (logStep <= 0) throw new IllegalArgumentException("deltaLog10 must be positive");
            this.logStep = logStep;
        }

        /**
         * Resample IM values onto a log-spaced grid, linearly interpolating
         * the corresponding response values.
         *
         * @param imValues original IM array (must be sorted ascending)
         * @param respValues original response array
         * @return transformed IM and response arrays
         */
        @Override
        public Result transform(double[] imValues, double[] respValues) {
            if (imValues.length != respValues.length) {
                throw new IllegalArgumentException("IM and response arrays must have same length");
            }
            double[] newIM = NumericUtil.distributeLogSpacedValues(
                    imValues[0],
                    imValues[imValues.length - 1],
                    logStep
            );

            double[] newResp = new double[newIM.length];

            for (int i = 0; i < newIM.length; i++) {
                newResp[i] = NumericUtil.linearInterpolate(
                        imValues,
                        respValues,
                        newIM[i]
                );
            }

            return new Result(newIM, newResp);
        }
    }
    

    /**
     * Transformer that resamples IM and response arrays onto a log-spaced IM grid.
     * Linear interpolation is used for the response values.
     */
    class CustomImTransformation implements ImValueTransformer {
    	
    	private final double [] customImValues;
    
    	
        public CustomImTransformation(double [] customImValues) {
            if (customImValues == null || customImValues.length == 0) {
                throw new IllegalArgumentException("customImValues cannot be null or empty");
            }

            for (int i = 0; i < customImValues.length; i++) {
                double v = customImValues[i];
                if (Double.isNaN(v) || Double.isInfinite(v)) {
                    throw new IllegalArgumentException(
                            "customImValues contains invalid number at index " + i);
                }
            }

            NumericUtil.checkMonotonicIM("customImValues", customImValues);

            this.customImValues = customImValues.clone();

        }
    	@Override
	    public Result transform(double[] imValues, double[] respValues) {
	
	        double[] newIM = customImValues.clone();
	        double[] newResp = new double[newIM.length];
	
	        for (int i = 0; i < newIM.length; i++) {
	            newResp[i] = NumericUtil.linearInterpolate(
	                    imValues,
	                    respValues,
	                    newIM[i]
	            );
	        }
	
	        return new Result(newIM, newResp);
	    }
	}

    /**
     * Placeholder transformer for future resampling logic.
     * Currently behaves as identity.
     */
    class PlaceholderImTransformation implements ImValueTransformer {
        @Override
        public Result transform(double[] imValues, double[] respValues) {
            // TODO implement optional interpolation (logspace, linspace, fixed points, etc.)
            return new Result(imValues.clone(), respValues.clone());
        }
    }
}
    
