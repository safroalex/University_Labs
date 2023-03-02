package lab1.strategies;

import lab1.common.Point;

public class Teleport implements MoveStrategy {
    @Override
    public Point getNewPos(Point start, Point dest) {
        return dest;
    }
}
