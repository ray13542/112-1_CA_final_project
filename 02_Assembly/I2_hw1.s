.data
    n: .word 10
    
.text
.globl __start

FUNCTION:
    # Todo: Define your own function in HW1
    # You should store the output into x10
    addi x11, x10, 0
    addi x10, x0, 0
T:  
    addi sp, sp, -8
    sw x11, 4(sp)
    sw ra, 0(sp)
    addi t1, x0, 1 #tmp1 = 1
    bgt x11, t1, else #iF(n <= 1)
    addi x10, x0, 2   #RESult = 2
    addi sp, sp, 8
    jalr x0, 0(ra) #back to else
else:
    srai x11, x11, 1 #n /= 2
    jal ra, T
    lw x11, 4(sp)   #n stored in a0
    lw ra, 0(sp)
    addi sp, sp, 8
    addi t2, x0, 5
    addi t3, x0, 6
    mul x10, x10, t2 #RESult = 5*RESult
    mul a6, t3, x11 #a6 = 6*n
    add x10, x10, a6 #RESult += 6n
    addi x10, x10, 4 #RESult += 4
    jalr x0, 0(ra)

# Do NOT modify this part!!!
__start:
    la   t0, n
    lw   x10, 0(t0)
    jal  x1,FUNCTION
    la   t0, n
    sw   x10, 4(t0)
    addi a0,x0,10
    ecall