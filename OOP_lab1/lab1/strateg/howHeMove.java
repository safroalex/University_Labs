package strateg;

import common.Point;
import java.security.InvalidParameterException;
import static java.lang.Math.sqrt;

public class howHeMove implements MoveStrategy {
    private final double avrVelocity;

    public howHeMove(double avrVelocity) throws InvalidParameterException { //protection from bad coders
        if (avrVelocity < 0) throw new InvalidParameterException("it must be more than zero");
        this.avrVelocity = avrVelocity;
    }

    @Override
    public Point getNewPos(Point start, Point stop) {
        double xIsLeft = stop.x - start.x;
        double yIsLeft = stop.y - start.y;
        double distance = sqrt(xIsLeft * xIsLeft + yIsLeft * yIsLeft);
        if (distance < avrVelocity) return stop;                              // we arrived
        return new Point (start.x + avrVelocity * xIsLeft / distance,      //next point
                start.y + avrVelocity * yIsLeft / distance);
    }
}
