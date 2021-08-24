from typing import Any, List


class SinglyLinkedList:

    def __init__(self):
        self._root = None
        self._tail = None
        self._count = 0

    def __len__(self):
        return self._count

    def append(self, item: Any) -> None:
        node = [item, None]

        if self._root is None:
            self._root = node
            self._tail = node
        else:
            self._tail[1] = node
            self._tail = node
        self._count += 1

    def extend(self, other) -> None:
        if other._root is None:
            return

        if self._root is None:
            self._root = other._root
            self._tail = other._tail
        else:
            self._tail[1] = other._root
            self._tail = other._tail
        self._count += len(other)

    def tolist(self) -> List[Any]:
        result = []
        curr = self._root
        while curr is not None:
            result.append(curr[0])
            curr = curr[1]
        return result
