from typing import Any


class SinglyLinkedList:

    def __init__(self):
        self._head = None
        self._tail = None
        self._count = 0

    def __len__(self):
        return self._count

    def __iter__(self):
        return SinglyLinkedListIterator(self)

    @property
    def head(self):
        return self._head

    def append(self, item: Any) -> None:
        node = [item, None]

        if self._head is None:
            self._head = node
            self._tail = node
        else:
            self._tail[1] = node
            self._tail = node
        self._count += 1

    def extend(self, other) -> None:
        assert isinstance(other, SinglyLinkedList)

        if other._head is None:
            return

        if self._head is None:
            self._head = other._head
            self._tail = other._tail
        else:
            self._tail[1] = other._head
            self._tail = other._tail
        self._count += len(other)


class SinglyLinkedListIterator:

    def __init__(self, linked_list):
        self._curr = linked_list.head

    def __iter__(self):
        return self

    def __next__(self):
        if self._curr is None:
            raise StopIteration
        else:
            item = self._curr[0]
            self._curr = self._curr[1]
            return item
