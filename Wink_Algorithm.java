import java.util.*;
//import org.json.*;

class Main {

    static double maxdegbearing = 60;
    static double maxmeterdist = 100;

    static double userlat = 0;
    static double userlon = 0;
    static double userheading = 120;



    public static void main(String[] args) {

        

        double temp_dist;
        double temp_bearingdiff;
        double[] temp_coords;
        double temp_cost;

        int best_index = 0;
        double best_cost = 3;

        double[][] data  = getRandomCoords();

        for (int i = 0; i < data.length; i++) {

            //temp_coords = getJSONCoordinates(i, JSONbuildings);
            temp_coords = data[i];
            System.out.println(Arrays.toString(temp_coords));
            
            temp_bearingdiff = calcBearingDifference(userheading, radToDegree(calcRadBearing(userlat, userlon, temp_coords[0],temp_coords[1])));

            if (temp_bearingdiff < maxdegbearing) {

                temp_dist = calcMeterDistance(userlat, userlon, temp_coords[0],temp_coords[1]);

                if(temp_dist < maxmeterdist) {

                    temp_cost = findCost(temp_dist, temp_bearingdiff);
                    System.out.println(temp_cost);

                    if (temp_cost < best_cost) {

                        best_index = i;
                        best_cost = temp_cost;
        
                    }
                }
            }   
        }

        System.out.println(best_index);
        System.out.println(best_cost);

    }

    public static double[][] getRandomCoords(){

        Random ra = new Random();
        double randomValue;

        double[][] randcoords = new double[150][2];

        for(int r = 0; r < 150; r++){

            for(int c = 0; c < 2; c++){

                randomValue = -1 + (1 + 1) * ra.nextDouble();
                randcoords[r][c] = randomValue;

                
            }         
        }
        return randcoords;

    }

    /* //returns double aray with index of 2, with the first index being latitude
    //and the second index being longitude
    public static double[] getJSONCoordinates(int index, JSONObject[] buildings){

        //assume we have an array of JSON objects representing individual
        //buildings on campus, as specified by NU Buildings API, called buildings.

    
        JSONObject req = new JSONObject(buildings[index]);
        double[] coords = new double[2];

        coords[0] = req.getDouble("lat");
        coords[1] = req.getDouble("lon");

        return coords;
    } */


    public static double findCost(double dist, double bearing) {

        return Math.pow(dist / maxmeterdist, 2) + Math.pow(bearing / maxdegbearing, 2);
    }

    private static double calcRadBearing(double lat1, double lon1, double lat2, double lon2) {

        double y = Math.sin(lon2 - lon1) * Math.cos(lat2);
        double x = Math.cos(lat1) * Math.sin(lat2) - Math.sin(lat1) * Math.cos(lat2) * Math.cos(lon2 - lon1);
        double radbearing = Math.atan2(y, x);

        if (radbearing < 0) {

            radbearing = 2 * Math.PI - radbearing;
        }


        return radbearing;


    }

    public static double calcMeterDistance(double lat1, double lon1, double lat2, double lon2) {
        double a = 6378137, b = 6356752.314245, f = 1 / 298.257223563; // WGS-84 ellipsoid params
        double L = Math.toRadians(lon2 - lon1);
        double U1 = Math.atan((1 - f) * Math.tan(Math.toRadians(lat1)));
        double U2 = Math.atan((1 - f) * Math.tan(Math.toRadians(lat2)));
        double sinU1 = Math.sin(U1), cosU1 = Math.cos(U1);
        double sinU2 = Math.sin(U2), cosU2 = Math.cos(U2);

        double sinLambda, cosLambda, sinSigma, cosSigma, sigma, sinAlpha, cosSqAlpha, cos2SigmaM;
        double lambda = L, lambdaP, iterLimit = 100;
        do {
            sinLambda = Math.sin(lambda);
            cosLambda = Math.cos(lambda);
            sinSigma = Math.sqrt((cosU2 * sinLambda) * (cosU2 * sinLambda) +
                (cosU1 * sinU2 - sinU1 * cosU2 * cosLambda) * (cosU1 * sinU2 - sinU1 * cosU2 * cosLambda));
            if (sinSigma == 0)
                return 0; // co-incident points
            cosSigma = sinU1 * sinU2 + cosU1 * cosU2 * cosLambda;
            sigma = Math.atan2(sinSigma, cosSigma);
            sinAlpha = cosU1 * cosU2 * sinLambda / sinSigma;
            cosSqAlpha = 1 - sinAlpha * sinAlpha;
            cos2SigmaM = cosSigma - 2 * sinU1 * sinU2 / cosSqAlpha;
            if (Double.isNaN(cos2SigmaM))
                cos2SigmaM = 0; // equatorial line: cosSqAlpha=0 (§6)
            double C = f / 16 * cosSqAlpha * (4 + f * (4 - 3 * cosSqAlpha));
            lambdaP = lambda;
            lambda = L + (1 - C) * f * sinAlpha *
                (sigma + C * sinSigma * (cos2SigmaM + C * cosSigma * (-1 + 2 * cos2SigmaM * cos2SigmaM)));
        } while (Math.abs(lambda - lambdaP) > 1e-12 && --iterLimit > 0);

        if (iterLimit == 0)
            return Double.NaN; // formula failed to converge

        double uSq = cosSqAlpha * (a * a - b * b) / (b * b);
        double A = 1 + uSq / 16384 * (4096 + uSq * (-768 + uSq * (320 - 175 * uSq)));
        double B = uSq / 1024 * (256 + uSq * (-128 + uSq * (74 - 47 * uSq)));
        double deltaSigma = B *
            sinSigma *
            (cos2SigmaM + B /
                4 *
                (cosSigma * (-1 + 2 * cos2SigmaM * cos2SigmaM) - B / 6 * cos2SigmaM *
                    (-3 + 4 * sinSigma * sinSigma) * (-3 + 4 * cos2SigmaM * cos2SigmaM)));
        double dist = b * A * (sigma - deltaSigma);

        return dist;
    }

    public static double radToDegree(double rad) {

        return rad * (180 / Math.PI);

    }

    public static double calcBearingDifference(double userheadingdeg, double degree){
        
        double degdifference = Math.abs(userheadingdeg-degree);

        return degdifference;
    }


}