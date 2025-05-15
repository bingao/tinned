// Generate all partitions of {1, ..., n} into exactly k non-empty subsets
fn generate_partitions(n: usize, k: usize) -> Vec<Vec<Vec<usize>>> {
    if n < k || k == 0 {
        return vec![]; // No valid partitions
    }

    if n == k {
        return vec![ (1..=n).map(|i| vec![i]).collect() ]; // Each element in its own group
    }

    if k == 1 {
        return vec![ vec![(1..=n).collect()] ]; // All elements in one group
    }

    let mut new_partitions = Vec::new();

    // Case 1: Use partitions from S(n-1, k), placing n into existing groups
    for partition in generate_partitions(n - 1, k) {
        for i in 0..partition.len() {
            let mut new_partition = partition.clone();
            new_partition[i].push(n);
            new_partitions.push(new_partition);
        }
    }

    // Case 2: Use partitions from S(n-1, k-1), creating a new singleton {n}
    for partition in generate_partitions(n - 1, k - 1) {
        let mut new_partition = partition.clone();
        new_partition.push(vec![n]);
        new_partitions.push(new_partition);
    }

    new_partitions
}

fn main() {
    let n = 5;
    let k = 2;
    let partitions = generate_partitions(n, k);

    println!("Partitions of {} into {} groups (size {}):", n, k, partitions.len());
    for p in &partitions {
        println!("{:?}", p);
    }
}
