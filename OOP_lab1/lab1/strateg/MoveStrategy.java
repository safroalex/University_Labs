package strateg;

import common.*;
public interface MoveStrategy {
    public Point getNewPos(Point start, Point dest);
}
