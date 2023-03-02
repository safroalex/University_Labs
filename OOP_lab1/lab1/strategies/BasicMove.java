package lab1.strategies;

import lab1.common.Point;

import java.security.InvalidParameterException;

import static java.lang.Math.*;

public class BasicMove implements MoveStrategy {
    private final double absVelocity;

    public BasicMove(double absVelocity) throws InvalidParameterException {
        if (absVelocity < 0) throw new InvalidParameterException("absVelocity should be above or equal to zero");
        this.absVelocity = absVelocity;
    }

    @Override
    public Point getNewPos(Point start, Point dest) {
        double xDist = dest.x - start.x;
        double yDist = dest.y - start.y;
        double distance = sqrt(xDist * xDist + yDist * yDist);
        if (distance < absVelocity) return dest;
        return new Point(start.x + absVelocity * xDist / distance, start.y + absVelocity * yDist / distance);
    }
}
