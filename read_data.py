def main():
    T10_count = 0
    T10_set = set()
    T40_count = 0
    T40_set = set()
    with open("dataset/T10I4D100K.txt") as file:
        data = file.readlines()
        for line in data:
            line = line.rstrip()
            for num in line.split(' '):
                T10_count += 1
                T10_set.add(num)

    with open("dataset/T40I10D100K.txt") as file:
        data = file.readlines()
        for line in data:
            line = line.rstrip()
            for num in line.split(' '):
                T40_count += 1
                T40_set.add(num)
    print("T10: ", T10_count, len(T10_set))
    print("T40: ", T40_count, len(T40_set))


if __name__ == '__main__':
    main()
