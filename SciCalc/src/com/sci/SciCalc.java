package com.sci;

/**
 * A scientific calculator library providing mathematical elements, operations, and constants
 * <p>
 * written by Angelos Hirmiz as a personal project
 * <p>
 * 
 * @author Angelos Hirmiz
 * @version 1.0
 */
public class SciCalc {

    // Prevent instantiation
    private SciCalc() {}

    /**
     * A collection of basic mathematical elements and operations
     */
    public static class Elem {

        // Prevent instantiation
        private Elem() {}
        
        /**
         * Accepts all real numbers and returns the magnitude of the input/the positive real value of the input
         * 
         * abs(num) = |num|
         * 
         * i.e., abs(5) = 5 : abs(-2) = 2 : abs(-0.707) = 0.707 : abs(0) = 0
         * 
         * @param num Any real number
         * @return The magnitude of num
         * 
         */
        public static double abs(double num) {
            return num < 0 ? -num : num; // Return -num if num < 0 or num if num > 0
        }

        /**
         * The exponential function, defined as the constant base e (2.7182818...) raised to a number
         * 
         * exp(num) = e^num
         * 
         * i.e., exp(1) = 2.7182818... : exp(0) = 1 : exp(3) = 20.0855369...
         * 
         * @param num Index of the base e
         * @return e raised to the index (e^index)
         */
        public static double exp(double num) {
            if (num == 0) return 1; // Return 1 for an index of 0, (anything raised to 0 is 1.0)

            else {
                final double PRECISION = 1E-15; // Threshold for which additional sums are negligible
                double sum = 1.0; // Value of the series expansion for n = 0
                double term = 1.0; // Value of initial term (n = 0)
                int n = 1; // Current term number

                while(abs(term) > PRECISION) { // Continuously sum terms until the term value is negligible
                    term *= num / n; // Update term value
                    sum += term; // Add term to overall sum
                    n++; // Increase term number
                }
                return sum;
            }
        }

        public static double log(double num) {
            if (num < 0) return Double.NaN;
            else if (num == 0) return Double.NEGATIVE_INFINITY;
            else {
                
            }
        }

        public static double floor(double num) {
            throw new UnsupportedOperationException("Not supported yet.");
        }
        
        /**
         * Defined as the method that raises a number to a power. 
         * 
         * pow(a, b) = a^b
         * 
         * Edge cases like 0^0 = 1, negative bases, etc are handled
         * 
         * i.e., pow(2, 2) = 4 : pow(9, 1/2) = 3 : pow(1.5, 1.5) = 1.837
         * 
         * @param base The base value.
         * @param index The index value.
         * @return The result of the base raised to the index, often expressed mathematically as
         * or NaN if the result is undefined numerically
         */
        public static double pow(double base, double index) {
            // Handle trivial cases
            if (index == 0) return 1.0; // Will indeed return 1 for 0^0 case by convention

            // Cases for base being 0, return 0 for index > 0 and undefined for negative indices
            if (base == 0) {
                if (index > 0) return 0.0;
                else return Double.NaN;
            }
            
            // Check for negative bases where the index is a non-integer
            if (base < 0 && index != floor(index)) return Double.NaN; 

            // Check for integer indices
            if (index == floor(index)) {
                double result = 1.0;
                int intIndex = (int) abs(index); // Ensure an integer index

                for (int i = 0; i < intIndex; i++) { // Compute by repeated multiplication
                    result *= base;
                }

                return index < 0 ? 1.0 / result : result; //Check for negative integer index, otherwise return the result
            }
            
            return exp(index * log(base)); // Use the identity base^index = exp(index * log(base))
        }

    }
    
    /**
     * A collection of scientific operations and functions
     */
    public static class Sci {
        // Prevent instantiation
        private Sci() {}
    }
    /**
     * A collection of important mathematical constants
     */
    public static class Const {
        // Prevent instantiation
        private Const() {}

        /**
         * The base of the natural logarithm, also known as Euler's number
         */
        public static final double E = 2.718281828459045;
        /**
         * The ratio of a circle's circumference to its diameter
         */
        public static final double PI = 3.141592653589793;
        /**
         * The golden ratio, often denoted by the Greek letter phi (φ)
         */
        public static final double PHI = 1.618033988749895;
        /**
         * The square root of 2, the length of the diagonal of a square with side length 1
         */
        public static final double SQRT2 = 1.4142135623730951;
        /**
         * The square root of 3, the length of the diagonal of a cube with side length 1
         */
        public static final double SQRT3 = 1.7320508075688772;
        /**
         * The natural logarithm of 2, the power to which e must be raised to obtain the value 2
         */
        public static final double LN2 = 0.6931471805599453;
        /**
         * The natural logarithm of 10, the power to which e must be raised to obtain the value 10
         */
        public static final double LN10 = 2.302585092994046;
        /**
         * The Euler-Mascheroni constant, which appears in various problems in number theory and analysis
         */
        public static final double EULER = 0.5772156649015329;
        
    }
}
