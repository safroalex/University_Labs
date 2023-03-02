package lab1.strategies;

import lab1.common.Point;

public interface MoveStrategy {
    public Point getNewPos(Point start, Point dest);
}
