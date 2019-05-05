var Main = (function () {
    function Main() {
    }
    Main.main = function (args) {
        var temp_dist;
        var temp_bearingdiff;
        var temp_coords;
        var temp_cost;
        var best_index = 0;
        var best_cost = 999999999.0;
        var data = Main.getRandomCoords();
        for (var i = 0; i < data.length; i++) {
            {
                temp_coords = data[i];
                console.info(temp_coords);
                temp_bearingdiff = Main.calcBearingDifference(Main.userheading, Main.radToDegree(Main.calcRadBearing(Main.userlat, Main.userlon, temp_coords[0], temp_coords[1])));
                //if (temp_bearingdiff < Main.maxdegbearing) {
                    temp_dist = Main.calcMeterDistance(Main.userlat, Main.userlon, temp_coords[0], temp_coords[1]);
                    //if (temp_dist < Main.maxmeterdist) {
                        temp_cost = Main.findCost(temp_dist, temp_bearingdiff);
                        console.info(temp_cost);
                        if (temp_cost < best_cost) {
                            best_index = i;
                            best_cost = temp_cost;
                        }
                    //}
                //}
            }
            ;
        }
        console.info(best_index);
        console.info(best_cost);
    };
    Main.getRandomCoords = function () {
        var randomValue;
        var randcoords = (function (dims) { var allocate = function (dims) { if (dims.length == 0) {
            return 0;
        }
        else {
            var array = [];
            for (var i = 0; i < dims[0]; i++) {
                array.push(allocate(dims.slice(1)));
            }
            return array;
        } }; return allocate(dims); })([150, 2]);
        for (var r = 0; r < 150; r++) {
            {
                for (var c = 0; c < 2; c++) {
                    {
                        randomValue = -1 + (1 + 1) * Math.random();
                        randcoords[r][c] = randomValue;
                    }
                    ;
                }
            }
            ;
        }
        return randcoords;
    };
    Main.findCost = function (dist, bearing) {
        return Math.pow(dist / Main.maxmeterdist, 2) + Math.pow(bearing / Main.maxdegbearing, 2);
    };
    /*private*/ Main.calcRadBearing = function (lat1, lon1, lat2, lon2) {
        var y = Math.sin(lon2 - lon1) * Math.cos(lat2);
        var x = Math.cos(lat1) * Math.sin(lat2) - Math.sin(lat1) * Math.cos(lat2) * Math.cos(lon2 - lon1);
        var radbearing = Math.atan2(y, x);
        if (radbearing < 0) {
            radbearing = 2 * Math.PI - radbearing;
        }
        return radbearing;
    };
    Main.calcMeterDistance = function (lat1, lon1, lat2, lon2) {
        var a = 6378137;
        var b = 6356752.314245;
        var f = 1 / 298.257223563;
        var L = (function (x) { return x * Math.PI / 180; })(lon2 - lon1);
        var U1 = Math.atan((1 - f) * Math.tan(/* toRadians */ (function (x) { return x * Math.PI / 180; })(lat1)));
        var U2 = Math.atan((1 - f) * Math.tan(/* toRadians */ (function (x) { return x * Math.PI / 180; })(lat2)));
        var sinU1 = Math.sin(U1);
        var cosU1 = Math.cos(U1);
        var sinU2 = Math.sin(U2);
        var cosU2 = Math.cos(U2);
        var sinLambda;
        var cosLambda;
        var sinSigma;
        var cosSigma;
        var sigma;
        var sinAlpha;
        var cosSqAlpha;
        var cos2SigmaM;
        var lambda = L;
        var lambdaP;
        var iterLimit = 100;
        do {
            {
                sinLambda = Math.sin(lambda);
                cosLambda = Math.cos(lambda);
                sinSigma = Math.sqrt((cosU2 * sinLambda) * (cosU2 * sinLambda) + (cosU1 * sinU2 - sinU1 * cosU2 * cosLambda) * (cosU1 * sinU2 - sinU1 * cosU2 * cosLambda));
                if (sinSigma === 0)
                    return 0;
                cosSigma = sinU1 * sinU2 + cosU1 * cosU2 * cosLambda;
                sigma = Math.atan2(sinSigma, cosSigma);
                sinAlpha = cosU1 * cosU2 * sinLambda / sinSigma;
                cosSqAlpha = 1 - sinAlpha * sinAlpha;
                cos2SigmaM = cosSigma - 2 * sinU1 * sinU2 / cosSqAlpha;
                if (isNaN(cos2SigmaM))
                    cos2SigmaM = 0;
                var C = f / 16 * cosSqAlpha * (4 + f * (4 - 3 * cosSqAlpha));
                lambdaP = lambda;
                lambda = L + (1 - C) * f * sinAlpha * (sigma + C * sinSigma * (cos2SigmaM + C * cosSigma * (-1 + 2 * cos2SigmaM * cos2SigmaM)));
            }
        } while ((Math.abs(lambda - lambdaP) > 1.0E-12 && --iterLimit > 0));
        if (iterLimit === 0)
            return NaN;
        var uSq = cosSqAlpha * (a * a - b * b) / (b * b);
        var A = 1 + uSq / 16384 * (4096 + uSq * (-768 + uSq * (320 - 175 * uSq)));
        var B = uSq / 1024 * (256 + uSq * (-128 + uSq * (74 - 47 * uSq)));
        var deltaSigma = B * sinSigma * (cos2SigmaM + B / 4 * (cosSigma * (-1 + 2 * cos2SigmaM * cos2SigmaM) - B / 6 * cos2SigmaM * (-3 + 4 * sinSigma * sinSigma) * (-3 + 4 * cos2SigmaM * cos2SigmaM)));
        var dist = b * A * (sigma - deltaSigma);
        return dist;
    };
    Main.radToDegree = function (rad) {
        return rad * (180 / Math.PI);
    };
    Main.calcBearingDifference = function (userheadingdeg, degree) {
        var degdifference = Math.abs(userheadingdeg - degree);
        return degdifference;
    };
    return Main;
}());
Main.maxdegbearing = 60;
Main.maxmeterdist = 100;
Main.userlat = 0;
Main.userlon = 0;
Main.userheading = 120;
Main["__class"] = "Main";
Main.main(null);
