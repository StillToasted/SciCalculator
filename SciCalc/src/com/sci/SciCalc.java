package com.sci;

/**
 * A scientific calculator library providing mathematical elements, operations, and constants.
 * <p>
 * written by Angelos Hirmiz as a personal project.
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

        //Common precision threshold for approximations and sums
        private final static double PRECISION = 1e-15;
        
        /**
         * Accepts all real numbers and returns the magnitude of the input/the positive real value of the input.
         * 
         * abs(x) = |x|
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
         * The square root of a value, defined as the value that when squared yields the original input.
         * This is calculated through Newton's method, also known as linear approximation, to approach
         * more and more accurate approximations of the square root.
         * 
         * sqrt(x) = √x
         * 
         * i.e., sqrt(4) = 2 : sqrt(2.25) = 1.5 : sqrt(2) = 1.414...
         * 
         * @param num A value
         * @return The square root of the value
         */
        public static double sqrt(double num) {
            // Handle trivial cases
            if (num < 0) return Double.NaN;
            if (num == 0) return 0;
            if (num == 1) return 1;
            
            long bits = Double.doubleToLongBits(num);    // Literally magic, I couldn't explain these two first lines if I tried.
            bits = (bits >> 1) + (0x1FF8000000000000L);  // I stole this from stack overflow, apparently it's a Quake 3 algorithm
            double root = Double.longBitsToDouble(bits); // for inverse square roots. The final line assigns an initial guess.
            
            double difference; // The value that modulates the guess
            do { // Loop continuously until the modulater is negligible
                difference = -((root * root) - num) / (2 * root);
                root += difference;
            } while(abs(difference) > PRECISION); 
            return root; // Return the square root                                                                         
        }

        /**
         * The cube root of a value, defined as the value that when cubed yields the original input.
         * This is calculated through Newton's method, also known as linear approximation, to approach
         * more and more accurate approximations of the cube root.
         * 
         * cbrt(x) = ∛x
         * 
         * i.e., cbrt(8) = 2 : cbrt(-216) = -6 : cbrt(0) = 0 : cbrt(42) = 3.476...
         * 
         * @param num A value
         * @return The cube root of the value
         */
        public static double cbrt(double num) {
            // Handle trivial cases
            if (num == 0) return 0;
            if (num == 1) return 1;
            if (num < 0) return -cbrt(-num);

            long bits = Double.doubleToLongBits(num);    // Again, literally magic, I don't know how this works or why, I also
            bits = bits / 3 + 0x2A9F7893E0000000L;       // stole this from stack overflow. The final line assigns an initial guess.
            double cbrt = Double.longBitsToDouble(bits);
            
            double difference; // The value that modulates the guess
            do { // Loop continuously until the modulater is negligible
                difference = -((cbrt * cbrt * cbrt) - num) / (3 * cbrt * cbrt);
                cbrt += difference;
            } while (abs(difference) > PRECISION);
            return cbrt;
        }

        /**
         * The exponential function, defined as the constant base e (2.7182818...) raised to a number. This function
         * is defined through the taylor series expansion of exp(x), where exp(x) = 1 + (x^2)/2! + (x^3)/3! + ...
         * 
         * exp(x) = e^x
         * 
         * i.e., exp(1) = 2.7182818... : exp(0) = 1 : exp(3) = 20.0855369...
         * 
         * @param num Index of the base e
         * @return e raised to the index (e^index)
         */
        public static double exp(double num) {
            if (num == 0) return 1; // Return 1 for an index of 0, (anything raised to 0 is 1.0)
            if (num < 0) return 1.0 / exp(-num); // Case for negative inputs

            double sum = 1.0; // Value of the series expansion for n = 0
            double term = 1.0; // Value of initial term (n = 0)
            int n = 1; // Current term number

            while(abs(term) > PRECISION) { // Continuously sum terms until the term value is negligible
                term *= num / n; // Update term value
                sum += term; // Add term to overall sum
                n++; // Increase term number
            } 
            return sum; // Final value of the exponential series
        }

        /**
         * The natural logarithm, described as the inverse function of the natural exponential function.
         * This logarithm has the index e, euler's number, valued at 2.7182818...
         * This is calculated through Newton's method, also known as linear approximation
         * to obtain more and more accurate approximations for the logarithm.
         * 
         * log(x) = ln(x)
         * 
         * i.e., log(2) = 0.6931... : log(e) = 1 : log(exp(10)) = 10
         * 
         * @param num The exponential yield, the value acquired when the logarithm is exponentiated
         * @return The index of the exponential function
         */
        public static double log(double num) {
            if (num < 0) return Double.NaN; // Case for negative inputs
            else if (num == 0) return Double.NEGATIVE_INFINITY; // General convention for log(0), approaching negative infinity
            else if (num == 1) return 0; // exp(0) = 1
            
            double logarithm = num > 1000 ? 10 : (num > 1 ? num / 2 : num); // Initial guess of logarithm
            double difference; // The value that modulates the guess
            do { // Loop continuously until the modulater is negligible
                difference = (num / exp(logarithm)) - 1;
                logarithm += difference;
            } while (abs(difference) > PRECISION);

        return logarithm; // Return the logarithm
        }

        /**
         * The baseless logarithm, described as the inverse function of an exponential function with an unspecified
         * base. This logarithm uses the change of base formula to compute the baseless logarithm, alongside the natural logarithm.
         * 
         * logb(base, num) = log_b(num)
         * 
         * i.e., logb(2, 0.75) = 1.6817 : logb(10, 1000) = 3 : logb(5, 25) = 2
         * 
         * @param base The desired logarithm base
         * @param num The exponential yield, the value acquired when the logarithm is exponentiated to the base
         * @return The index of the exponential
         */
        public static double logb(double base, double num) {
            // Log bases cannot be 1, or <= 0. The input must be greater than or equal to 0
            if (base <= 0 || base == 1 || num < 0) return Double.NaN;
            if (num == 1) return 0.0; // Avoid redundant computation
            return log(num) / log(base); // Change of base formula
        }

        /**
         * In radians
         * 
         * The sine function, described as the y-coordinate of any point on the unit circle subtended
         * by an angle in radians from the positive x-axis. This is computed through the taylor series
         * expansion of sin(x). The sine function has the property that sin(-x) = -sin(x)
         * 
         * i.e., sin(2*PI) = 0 : sin(PI/6) = 0.5 : sin(-PI / 4) = -0.7071...
         * 
         * @param angle An angle (radians) from the positive x-axis on the unit-circle
         * @return The y-value of the point on the unit circle relative to the angle
         */
        public static double sin(double angle) {
            if (angle < 0) return -sin(-angle); // Sine is an odd function, return -sin(x) for sin(-x)
            if (angle > 2 * Const.PI) angle %= 2 * Const.PI; // Limit the angle 
            // Trivial cases
            if (angle == 0 || angle == 2 * Const.PI || angle == Const.PI) return 0;
            if (angle == Const.PI / 2) return 1;
            if (angle == (3 / 2) * Const.PI) return -1;
            
            double term = angle; // Current term of the series expansion
            double sum = term; // Current sum
            int n = 1; //Current term number

            while(abs(term) > PRECISION) { // Loop continously until the modulator is negligible
                term *= (angle * angle * -1) / ((2 * n + 1) * (2 * n)); // Update the term
                sum += term; // Update the sum
                n++; // Increment the term number
            }
            return sum; // Final value of the sine series
        }

        /**
         * In radians or degrees
         * 
         * The sine function, described as the y-coordinate of any point on the unit circle subtended
         * by an angle in either radians or degrees. This is computed through the taylor series expansion of sin(x). The sine function
         * has the property that sin(-x) = -sin(x)
         * 
         * i.e., sin(2*PI) = 0 : sin(PI/6) = 0.5 : sin(-PI / 4) = -0.7071...
         * 
         * or
         * 
         * i.e., sin(360 deg) = 0 : sin(30 deg) = 0.5 : sin(-45 deg) = -0.7071...
         * 
         * @param angle An angle (radians/degrees) from the positive x-axis on the unit-circle
         * @return The y-value of the point on the unit circle relative to the angle
         */
        public static double sin(double angle, boolean inDegrees) {
            if (inDegrees) return sin((Const.PI / 180) * angle); // Accept an input in degrees and convert back to radians
            else return sin(angle); // Return the regular function if !inDegrees
        }

        /**
         * In radians
         * 
         * The cosine function, described as the x-coordinate of any point on the unit circle subtended
         * by an angle in radians from the positive x-axis. This is computed through the taylor series
         * expansion of cos(x). The cosine function has the property that cos(-x) = cos(x)
         * 
         * i.e., cos(PI/2) = 0 : cos(PI) = -1 : cos(PI/6) = 0.866...
         * 
         * @param angle An angle (radians) from the positive x-axis on the unit-circle
         * @return The x-value of the point on the unit circle relative to the angle
         */
        public static double cos(double angle) {
            if (angle < 0) return cos(-angle); // Cosine is an even function, return cos(x) for cos(-x)
            if (angle > 2 * Const.PI) angle %= 2 * Const.PI; // Limit the angle 
            // Trivial cases
            if (angle ==  Const.PI / 2 || angle == Const.PI * (3.0 / 2.0)) return 0;
            if (angle == 0) return 1;
            if (angle == Const.PI) return -1;

            double term = 1; // Current term of the series expansion
            double sum = term; // Current sum
            int n = 1; //Current term number

            while(abs(term) > PRECISION) { // Loop continuously until loop modulator is negligible
                term *= (angle * angle * -1) / ((2 * n - 1) * (2 * n)); // Update the term
                sum += term; // Update the sum
                n++; // Increment the term number
            }
            return sum; // Final value of the cosine series
        }

        /**
         * In radians or degrees
         * 
         * The cosine function, described as the x-coordinate of any point on the unit circle subtended
         * by an angle in radians or degrees from the positive x-axis. This is computed through the taylor series
         * expansion of cos(x). The cosine function has the property that cos(-x) = cos(x)
         * 
         * i.e., cos(PI/2) = 0 : cos(PI) = -1 : cos(PI/6) = 0.866...
         * 
         * or
         * 
         * i.e., cos(90 deg) = 0 : cos(180 deg) = -1 : cos(30 deg) = 0.866...
         * 
         * @param angle An angle (radians/degrees) from the positive x-axis on the unit-circle
         * @return The x-value of the point on the unit circle relative to the angle
         */
        public static double cos(double angle, boolean isDegrees) {
            if (isDegrees) return cos((Const.PI / 180) * angle); // Accept an input in degrees and convert back to radians
            else return cos(angle); // Return the regular function if !inDegrees
        }

        /**
         * The tangent function, described as the ratio of the y and x coordinates of any point on the unit circle
         * subtended by an angle in radians from the positive x-axis. This is computed through the identity that
         * tan(x) = sin(x) / cos(x), where sine and cosine are evaluated through their taylor series expansions. The
         * tangent function has the property that tan(-x) = -tan(x)
         * 
         * i.e., tan(PI) = 0 : tan(PI/4) = 1 : tan(-PI/8) = -0.4142...
         * 
         * @param angle An angle (radians) from the positive x-axis on the unit-circle
         * @return The value of y / x, where x and y are the coordinates of a point on the unit circle relative to the angle
         */
        public static double tan(double angle) {
            if (angle < 0) return -tan(-angle); // Tangent is an odd function, return -tan(x) for tan(-x)
            if (angle > 2 * Const.PI) angle %= 2 * Const.PI; // Limit the angle
            if (cos(angle) == 0) return Double.NaN; // Case when the identity sin(x) / cos(x) becomes undefined

            return sin(angle) / cos(angle); // Traditional identity for tan(x) = sin(x) / cos(x)
        }

        /**
         * The tangent function, described as the ratio of the y and x coordinates of any point on the unit circle
         * subtended by an angle in radians or degrees from the positive x-axis. This is computed through the identity that
         * tan(x) = sin(x) / cos(x), where sine and cosine are evaluated through their taylor series expansions. The
         * tangent function has the property that tan(-x) = -tan(x)
         * 
         * i.e., tan(PI) = 0 : tan(PI/4) = 1 : tan(-PI/8) = -0.4142...
         * 
         * or
         * 
         * i.e., tan(180 deg) = 0 : tan(45 deg) = 1 : tan(-22.5 deg) = -0.4142...
         * 
         * @param angle An angle (radians/degrees) from the positive x-axis on the unit-circle
         * @return The value of y / x, where x and y are the coordinates of a point on the unit circle relative to the angle
         */
        public static double tan(double angle, boolean isDegrees) {
            if (isDegrees) return tan((Const.PI / 180) * angle); // Accept an input in degrees and convert back to radians
            else return (tan(angle)); // Return the regular function if !inDegrees
        }

        /**
         * The cosecant function, traditionally described as the reciprocal of the sine function, or 1.0 / sin(x). Is undefined
         * when sin(x) = 0
         * 
         * i.e., csc(PI/4) = 1.414... : csc(-PI/6) = -2 : csc(PI/2) = 1
         * 
         * @param angle An angle (radians) from the positive x-axis on the unit-circle
         * @return The value of 1 / y, where y is the second coordinate of any point on the unit-circle
         */
        public static double csc(double angle) {
            if (sin(angle) == 0) return Double.NaN; // Case when 1 / sin(x) is undefined
            else return 1.0 / sin(angle); // Use identity that csc(x) = 1 / sin(x)
        }
        
        /**
         * The cosecant function, traditionally described as the reciprocal of the sine function, or 1.0 / sin(x). Is undefined
         * when sin(x) = 0
         * 
         * i.e., csc(PI/4) = 1.414... : csc(-PI/6) = -2 : csc(PI/2) = 1
         * 
         * or
         * 
         * i.e., csc(45 deg) = 1.414... : csc(-30 deg) = -2 : csc(90 deg) = 1 
         * 
         * @param angle An angle (radians/degrees) from the positive x-axis on the unit-circle
         * @return The value of 1 / y, where y is the second coordinate of any point on the unit-circle
         */
        public static double csc(double angle, boolean isDegrees) {
            if (isDegrees) return csc((Const.PI / 180) * angle); // Accept an input in degrees and convert back to radians
            else return csc(angle); // Return the regular function if !inDegrees
        }

        /**
         * The secant function, traditionally described as the reciprocal of the cosine function, or 1.0 / cos(x). Is undefined
         * when cos(x) = 0
         * 
         * i.e., sec(PI/4) = 1.414... : sec(-PI/3) = 2 : sec(PI) = -1
         * 
         * @param angle An angle (radians) from the positive x-axis on the unit-circle
         * @return The value of 1 / x, where x is the second coordinate of any point on the unit-circle
         */
        public static double sec(double angle) {
            if (cos(angle) == 0) return Double.NaN; // Case when 1 / cos(x) is undefined
            else return 1.0 / cos(angle); // Use identity that sec(x) = 1 / cos(x)
        }

        /**
         * The secant function, traditionally described as the reciprocal of the cosine function, or 1.0 / cos(x). Is undefined
         * when cos(x) = 0
         * 
         * i.e., sec(PI/4) = 1.414... : sec(-PI/3) = 2 : sec(PI) = -1
         * 
         * or
         * 
         * i.e., sec(45 deg) = 1.414... : sec(-60 deg) = 2 : sec(180 deg) = -1
         * 
         * @param angle An angle (radians/degrees) from the positive x-axis on the unit-circle
         * @return The value of 1 / x, where x is the second coordinate of any point on the unit-circle
         */
        public static double sec(double angle, boolean isDegrees) {
            if (isDegrees) return sec((Const.PI / 180) * angle); // Accept an input in degrees and convert back to radians
            else return sec(angle); // Return the regular function if !inDegrees
        }

        /**
         * The cotangent function, traditionally described as the reciprocal of the tangent function, or 1.0 / tan(x). Is undefined
         * when sin(x) = 0
         * 
         * i.e., cot(PI/4) = 1 : cot(-PI/8) = -2.414... : cot(PI/2) = 0
         * 
         * @param angle An angle (radians) from the positive x-axis on the unit-circle
         * @return The value of x / y, where x and y are the coordinates of a point on the unit circle relative to the angle
         */
        public static double cot(double angle) {
            if (tan(angle) == 0) return Double.NaN; // Case when 1 / cos(x) is undefined
            else return 1.0 / tan(angle); // Use identity that sec(x) = 1 / cos(x)
        }

        /**
         * The cotangent function, traditionally described as the reciprocal of the tangent function, or 1.0 / tan(x). Is undefined
         * when sin(x) = 0
         * 
         * i.e., cot(PI/4) = 1 : cot(-PI/8) = -2.414... : cot(PI/2) = 0
         * 
         * or
         * 
         * i.e., cot(45 deg) = 1 : cot(-22.5 deg) = -2.414... : cot(90 deg) = 0
         * 
         * @param angle An angle (radians/degrees) from the positive x-axis on the unit-circle
         * @return The value of x / y, where x and y are the coordinates of a point on the unit circle relative to the angle
         */
        public static double cot(double angle, boolean isDegrees) {
            if (isDegrees) return cot((Const.PI / 180) * angle); // Accept an input in degrees and convert back to radians
            else return cot(angle); // Return the regular function if !inDegrees
        }

        /**
         * The arcsine function, inverse to the ordinary sine function, returning an angle in radians. Has the property:
         * arcsin(-x) = -arcsin(x); is restricted by a domain between -1 and 1 and range between -PI/2 and PI/2
         * 
         * for sin(x) = a, a = asin(x)
         * 
         * i.e., asin(1) = 1.57... : asin(0.707) = 0.785... : asin(-0.5) = -0.523... : asin(0) = 0
         * 
         * @param len The length of the terminal arm on the unit circle
         * @return The angle that is subtended between the positive x-axis and the terminal arm in radians
         */
        public static double asin(double len) {
            if (len < 0) return -asin(-len); // Arcsine is an odd function, return -asin(x) for asin(-x)
            // Trivial cases
            if (len == 0) return 0;
            if (len > 1) return Double.NaN;
            if (len == 1) return Const.PI / 2;

            double angle = (Const.PI / 2) - sqrt(2 * (1 - len)); // Series based initial guess of arcsine
            double difference;
            do { // Loop continuously until loop modulator is negligible
                difference = -(sin(angle) - len) / (cos(angle));
                angle += difference;
            } while(abs(difference) > PRECISION);
            return angle; // Return the final angle
        }

        /**
         * The arcsine function, inverse to the ordinary sine function, returning an angle in radians/degrees. Has the 
         * property: arcsin(-x) = -arcsin(x); is restricted by a domain between -1 and 1 and range between -PI/2 and PI/2 or
         * -90 and 90
         * 
         * for sin(x) = a, a = asin(x)
         * 
         * i.e., asin(1) = 1.57... : asin(0.707) = 0.785... : asin(-0.5) = -0.523... : asin(0) = 0
         * 
         * or
         * 
         * i.e., asin(1) = 90 : asin(0.707) = 45 : asin(-0.5) = -30 : asin(0) = 0
         * 
         * @param len The length of the terminal arm on the unit circle
         * @return The angle that is subtended between the positive x-axis and the terminal arm in radians/degrees
         */
        public static double asin(double len, boolean isDegrees) {
            if (isDegrees) return asin(len) * 180 / Const.PI; // Convert the radian result to degrees
            else return asin(len); // If !degrees, return the regular function in radians
        }
        /**
         * The arccosine function, inverse to the ordinary cosine function, returning an angle in radians.
         * Is restricted by a domain between -1 and 1 and range between PI and 0 or
         * 180 and 0
         * 
         * for cos(x) = a, a = acos(x)
         * 
         * i.e., acos(1) = 0 : arcsin(0.707) = 0.785... : acos(-0.5) = 2.0943... : acos(0) = 1.57...
         * 
         * @param len The length of the terminal arm on the unit circle
         * @return The angle that is subtended between the positive x-axis and the terminal arm in radians/degrees
         */
        public static double acos(double len) {
            if (len == -1) return Const.PI;
            if (len == 0) return Const.PI / 2;
            if (len == 1) return 0;

            double angle = (len >= 0) ? 
            (Const.PI / 2) - sqrt(2 * (1 - len)) : // For len in [0, 1]
            Const.PI - sqrt(2 * (1 + len)); // For len in [-1, 0)
            double difference;
            do { // Loop continuously until loop modulator is negligible
                difference = (cos(angle) - len) / (sin(angle));
                angle += difference;
            } while (abs(difference) > PRECISION);
            return angle; // Return the final angle
        }
        
        /**
         * The arccosine function, inverse to the ordinary cosine function, returning an angle in radians.
         * Is restricted by a domain between -1 and 1 and range between PI and 0 or
         * 180 and 0
         * 
         * for cos(x) = a, a = acos(x)
         * 
         * i.e., acos(1) = 0 : arcsin(0.707) = 0.785... : acos(-0.5) = 2.0943... : acos(0) = 1.57...
         * 
         * @param len The length of the terminal arm on the unit circle
         * @return The angle that is subtended between the positive x-axis and the terminal arm in radians/degrees
         */
        public static double acos(double len, boolean isDegrees) {
            if (isDegrees) return acos(len) * 180 / Const.PI; // Convert the radian result to degrees
            else return acos(len); // If !degrees, return the regular function in radians
        }

        /**
         * The floor function, taking a number argument and yielding the lesser of the two integers that it falls between.
         * 
         * floor(x) = ⌊x⌋
         * 
         * i.e., floor(2.8) = 2 : floor(0.999) = 0 : floor(-10) = -10 : floor(-0.00001) = -1
         * 
         * @param num A value
         * @return The lesser of the two integers that the value falls between
         */
        public static double floor(double num) {
            int intPart = (int) num; // Obtain the integer component

            // If the number is already an integer, or positive, return the truncated value. Otherwise subtract 1 from the integer component
            return (num == intPart || num >= 0) ? intPart : intPart - 1; 
        }

        /**
         * The ceiling function, taking a number argument and yielding the greater of the two integers that it falls between.
         * 
         * ceil(x) = ⌈x⌉
         * 
         * i.e., ceil(2.8) = 3 : ceil(0.0001) = 1 : ceil(-10) = -10 : ceil(-0.5) = 0
         * 
         * @param num A value
         * @return The greater of the two integers that the value falls between
         */
        public static double ceil(double num) {
            int intPart = (int) num; // Obtain the integer component

            // If the number is already an integer, or negative, return the truncated value. Otherwise add 1 from the integer component
            return (num == intPart || num <= 0) ? intPart : intPart + 1;
        }

        /**
         * The standard rounding function, accepting all real numbers and returning the number rounded to the nearest integer. If the decimal component of
         * the argument is greater than 0.5, the argument is rounded up, and vice versa.
         * 
         * i.e., roundInt(2) = 2 : roundInt(-1.9) = -2 : roundInt(4.5) = 5
         * 
         * @param num A value
         * @return The value rounded to the nearest integer
         */
        public static double roundInt(double num) {
            if (num < 0) return -roundInt(-num); // Round numbers symmetrically

            int intPart = (int) num; // Obtain the integer component of the number
            double frac = num - intPart; // Obtain the decimal component of the number

            return (frac < 0.5) ? intPart : intPart + 1; // Round the number
        }

        /**
         * The extended rounding function, accepting all real numbers and a specified number of decimal places, for the number to be rounded to. If the decimal
         * place's value exceeds 5, the argument is rounded up to the nearest specified decimal place, and vice versa.
         * 
         * i.e., round(2.52, 0) = 3 : round(-4.55, 1) = -4.6 : round(-0.469, 2) = -0.47
         * 
         * @param num A value
         * @param decimalPlaces The number of decimal places the value is to be rounded
         * @return The value rounded to the specified number of decimal places
         */
        public static double round(double num, int decimalPlaces) {
            double factor = pow(10, decimalPlaces); // Compute the multiplication factor
            return roundInt(num * factor) / factor; // Round the number multiplied by the factor, then divide by the factor to obtain the rounded number
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
            if (base < 0 && index != (int) index) return Double.NaN; 

            // Check for integer indices
            if (index == (int) index) {
                double result = 1.0;
                int intIndex = (int) abs(index); // Ensure a positive integer index

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
    public static class Ops {
        // Prevent instantiation
        private Ops() {}
    }
    /**
     * A collection of important mathematical constants
     */
    public static class Const {
        // Prevent instantiation
        private Const() {}

        /**
         * The base of the natural logarithm or exponential function, also known as Euler's number
         */
        public static final double E = 2.718281828459045;
        /**
         * The ratio of a circle's circumference to its diameter
         */
        public static final double PI = 3.141592653589793;
        /**
         * The alternative circle constant, the ratio fo a circle's circumference to
         * its radius
         */
        public static final double TAU = 6.283185307179586;
        /**
         * The golden ratio, often denoted by the Greek letter phi
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
