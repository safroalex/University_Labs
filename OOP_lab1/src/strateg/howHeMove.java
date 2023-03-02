package strateg;

import startHero.Point;

import java.security.InvalidParameterException;

import static java.lang.Math.sqrt;

public class howHeMove implements start_stop_point {
    private final double avrgVelocity;

    public howHeMove(double avrgVelocity) throws InvalidParameterException {
        if (avrgVelocity < 0) throw new InvalidParameterException("it must be more than zero");
        this.avrgVelocity = avrgVelocity;
    }

    @Override
    public Point newPos(Point start, Point stop) {
        double xIsLeft = stop.x - start.x;
        double yIsLeft = stop.y - start.y;
        double distance = sqrt(xIsLeft * xIsLeft + yIsLeft * yIsLeft);
        if (distance < avrgVelocity) return stop;                           // we arrived
        return new Point (start.x + avrgVelocity * xIsLeft / distance,      //next point
                start.y + avrgVelocity * yIsLeft / distance);
    }
}
