package lab1.strategies;

import lab1.common.Point;

public class Stay implements MoveStrategy {
    @Override
    public Point getNewPos(Point start, Point dest) {
        return start;
    }
}
